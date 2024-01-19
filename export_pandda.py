# Copyright (c) 2024, Tobias Krojer, MAX IV Laboratory
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import pandas as pd
import gemmi

import getopt
import glob
import sys
import os
from datetime import datetime

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'lib')))
import processlib

def read_pandda_analyse_events_as_df(logger, panddaDir):
    eventdf = None
    if os.path.isfile(os.path.join(panddaDir, 'results', 'pandda_inspect_events.csv')):
        eventcsv = os.path.join(panddaDir, 'results', 'pandda_inspect_events.csv')
        logger.info('reading {0!s} as dataframe...'.format(eventcsv))
        eventdf = pd.read_csv(eventcsv)
    elif os.path.isfile(os.path.join(panddaDir, 'analyses', 'pandda_inspect_events.csv')):
        eventcsv = os.path.join(panddaDir, 'analyses', 'pandda_inspect_events.csv')
        logger.info('reading {0!s} as dataframe...'.format(eventcsv))
        eventdf = pd.read_csv(eventcsv)
    else:
        logger.error('cannot find {0!s}'.format(eventcsv))
    return eventdf

def get_refinement_method(logger,confidence, d, sample, event):
    if d[sample] == True:
        if "2fofc" in confidence:
            logger.info('{0!s} - {0!s}: at least one ligand clear in 2fofc map; setting ensemble refinement to False'.format(sample, event))
            d[sample] = False
        else:
            logger.info('{0!s} - {0!s}: ligand only visible in event map; leaving ensemble refinement as True'.format(sample, event))
    else:
        if "2fofc" in confidence:
            logger.info('{0!s} - {0!s}: ligand clear in 2fofc map; setting ensemble refinement to False'.format(sample, event))
            d[sample] = False
        else:
            logger.info('{0!s} - {0!s}: ligand only visible in event map; setting ensemble refinement as True'.format(sample, event))
            d[sample] = True
    return d

def get_sample_dict_for_export(logger, panddaDir):
    eventdf = read_pandda_analyse_events_as_df(logger, panddaDir)
    d = {}
    for index, row in eventdf.iterrows():
        sample = row['dtag']
        event = row['event_num']
        confidence = row['Ligand Confidence']
        placed = row['Ligand Placed']
        if placed:
            if sample not in d:
                d[sample] = []
            d = get_refinement_method(logger,confidence, d, sample, event)
    return d

def prepare_sample_refine_folder(logger, fragmaxDir, sample):
    os.chdir(fragmaxDir)
    if not os.path.isdir('5-refine'):
        logger.info("creating '5-refine' folder in {0!s}".format(fragmaxDir))
        os.mkdir('5-refine')
    os.chdir('5-refine')
    if not os.path.isdir(sample):
        logger.info("creating '{0!s}' folder in {1!s}/5-refine".format(sample, fragmaxDir))
        os.mkdir(sample)
    else:
        logger.info("'{0!s}' folder exists in {1!s}/5-refine".format(sample, fragmaxDir))

def sample_refine_folder_is_empty(logger, fragmaxDir, sample):
    sample_folder_empty = True
    os.chdir(os.path.join(fragmaxDir, '5-refine', sample))
    for f in glob.glob('*'):
        logger.warning("'{0!s}' folder in {1!s}/5-refine in not empty".format(sample, fragmaxDir))
        sample_folder_empty = False
    if sample_folder_empty:
        logger.info("'{0!s}' folder in {1!s}/5-refine in empty".format(sample, fragmaxDir))
    return sample_folder_empty


def backup_existing_files_folders_in_sample_dir(logger, fragmaxDir, sample):
    now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    logger.warning("'{0!s}' folder in {1!s}/5-refine in not empty; backing up files/ folders in 'backup_{2!s}'".format(sample, fragmaxDir, now))
    os.chdir(os.path.join(fragmaxDir, '5-refine', sample))
    os.mkdir('backup_' + now)
    for f in glob.glob('*'):
        if not f.startswith('backup_'):
            os.system('/bin/mv {0!s} backup_{1!s}'.format(f, now))

def prepare_ensemble_model(logger, panddaDir, sample):
    logger.info('{0!s} - preparing ensemble model'.format(sample))
    os.chdir(os.path.join(panddaDir, 'processed_datasets', sample))
    logger.info('changing directory to {0!s}'.format(os.path.join(panddaDir, 'processed_datasets', sample)))
    if os.path.isfile('merge_conformations.pdb'):
        logger.warning("removing 'merge_conformations.pdb'")
        os.system('/bin/rm merge_conformations.pdb')
    input_model = os.path.join(panddaDir, 'processed_datasets', sample, sample + '-pandda-input.pdb')
    if os.path.isfile(input_model):
        cmd = (
        'module load gopresto CCP4\n'
        'giant.merge_conformations modelled_structures/{0!s}-pandda-model.pdb {1!s}-pandda-input.pdb > /dev/null\n'.format(sample, sample)
           )
        logger.info("running 'giant.merge_conformations modelled_structures/{0!s}-pandda-model.pdb {1!s}-pandda-input.pdb'".format(sample, sample))
        os.system(cmd)
    else:
        logger.error('input model for pandda does not exist: {0!s}'.format(input_model))

def set_ensemble_refinement_mark(logger, fragmaxDir, sample, ensemble):
    logger.info("creating empty file 'REFINE_AS_ENSEMBLE' to indicate ensemble refinement with giant.quick_refine")
    os.chdir(os.path.join(fragmaxDir, '5-refine', sample))
    os.system('touch REFINE_AS_ENSEMBLE')

def copy_process_files(logger, fragmaxDir, sample):
    logger.info('copying process files from {0!s}'.format(os.path.join(fragmaxDir, '1-process', sample)))
    os.chdir(os.path.join(fragmaxDir, '5-refine', sample))
    os.system('/bin/cp {0!s}/process.* .'.format(os.path.join(fragmaxDir, '1-process', sample)))
    os.system('/bin/cp {0!s}/info.json .'.format(os.path.join(fragmaxDir, '1-process', sample)))

def copy_ligand_restraints(logger, fragmaxDir, sample, panddaDir):
    logger.info('copying ligand_files from {0!s}'.format(os.path.join(panddaDir, 'processed_datasets', sample)))
    os.chdir(os.path.join(fragmaxDir, '5-refine', sample))
    os.system('/bin/cp -R {0!s}/ligand_files .'.format(os.path.join(panddaDir, 'processed_datasets', sample)))

def copy_free_mtz(logger, fragmaxDir, sample):
    logger.info('copying free.mtz from {0!s}'.format(os.path.join(fragmaxDir, '2-initial_refine', sample)))
    os.chdir(os.path.join(fragmaxDir, '5-refine', sample))
    os.system('/bin/cp {0!s}/free.mtz .'.format(os.path.join(fragmaxDir, '2-initial_refine', sample)))
    os.system('ln -s init.mtz refine.mtz')

def copy_init_refine_files(logger, fragmaxDir, sample):
    logger.info('copying initial refinement files from {0!s}'.format(os.path.join(fragmaxDir, '2-initial_refine', sample)))
    os.chdir(os.path.join(fragmaxDir, '5-refine', sample))
    os.system('/bin/cp {0!s}/init.* .'.format(os.path.join(fragmaxDir, '2-initial_refine', sample)))

def copy_pdb_file(logger, fragmaxDir, sample, panddaDir, ensemble):
    os.chdir(os.path.join(fragmaxDir, '5-refine', sample))
    if ensemble:
        pdb = os.path.join(panddaDir, 'processed_datasets', sample, 'merge_conformations.pdb')
        par = os.path.join(panddaDir, 'processed_datasets', sample, 'merge_conformations-restraints.refmac.params')
        logger.info('copying ensemble PDB file {0!s}'.format(pdb))
        os.system('/bin/cp {0!s} .'.format(pdb))
        logger.info('copying refmac restraints PDB file {0!s}'.format(par))
        os.system('/bin/cp {0!s} .'.format(par))
        os.system('ln -s merge_conformations.pdb refine.pdb')
    else:
        pdb = os.path.join(panddaDir, 'processed_datasets', sample, 'modelled_structures', sample + '-pandda-model.pdb')
        logger.info('copying single conformer PDB file {0!s}'.format(pdb))
        os.system('/bin/cp {0!s} .'.format(pdb))
        os.system('ln -s {0!s}-pandda-model.pdb refine.pdb'.format(sample))

def assign_modelled_ligand_to_event_coordinate(logger, sample, panddaDir, eventdf):
    model = os.path.join(panddaDir, 'processed_datasets', sample, 'modelled_structures', sample + '-pandda-model.pdb')
    logger.info('trying to assign event (maps) to modelled ligands of type LIG in {0!s}'.format(model))
    structure = gemmi.read_structure(model, merge_chain_parts=True)
    matching_event_maps = []
    for mod in structure:
        for chain in mod:
            for residue in chain:
                if residue.name == 'LIG':
                    c = gemmi.Chain(chain.name)
                    c.add_residue(residue, 0)
                    for index, row in eventdf.iterrows():
                        dtag = row['dtag']
                        if dtag == sample:
                            eventidx = row['event_num']
                            x = row['x']
                            y = row['y']
                            z = row['z']
                            event = gemmi.Position(float(x), float(y), float(z))
                            distance = round(event.dist(c.calculate_center_of_mass()), 2)
                            lig = residue.name + '-' + chain.name + '-' + str(residue.seqid.num)
                            logger.info('event: {0!s} - x={1!s}, y={2!s}, z={3!s} - ligand: {4!s} - distance: {5!s}'.format(eventidx, x, y, z, lig, distance))
                            if distance < 8:
                                os.chdir(os.path.join(panddaDir, 'processed_datasets', sample))
                                for emap in glob.glob('*-pandda-output-event-*.mtz'):
                                    n_map = emap.replace('.mtz', '').split('-')[len(emap.replace('.mtz', '').split('-'))-1]
                                    if int(eventidx) == int(n_map):
                                        logger.info('matching event map is {0!s}'.format(emap))
                                        matching_event_maps.append(emap)
                                logger.info('ligand is within 8A of event, assuming they match...')
                            else:
                                logger.warning('ligand is more than 8A away from event; ignoring it...')
    if matching_event_maps:
        logger.info('found at least one matching event map for modelled ligand(s)')
    else:
        logger.error('could not find any matching event maps...')
    return matching_event_maps

def copy_event_maps(logger, fragmaxDir, sample, panddaDir, eventdf):
    matching_event_maps = assign_modelled_ligand_to_event_coordinate(logger, sample, panddaDir, eventdf)
    os.chdir(os.path.join(fragmaxDir, '5-refine', sample))
    for emap in matching_event_maps:
        emtz = os.path.join(panddaDir, 'processed_datasets', sample, emap)
        logger.info('copying event mtz file {0!s}'.format(emtz))
        os.system('/bin/cp {0!s} .'.format(emtz))

def copy_files(logger, panddaDir, sample, fragmaxDir, ensemble):
    eventdf = read_pandda_analyse_events_as_df(logger, panddaDir)
    if ensemble:
        prepare_ensemble_model(logger, panddaDir, sample)
        set_ensemble_refinement_mark(logger, fragmaxDir, sample, ensemble)
    copy_process_files(logger, fragmaxDir, sample)
    copy_ligand_restraints(logger, fragmaxDir, sample, panddaDir)
    copy_free_mtz(logger, fragmaxDir, sample)
    copy_init_refine_files(logger, fragmaxDir, sample)
    copy_pdb_file(logger, fragmaxDir, sample, panddaDir, ensemble)
    copy_event_maps(logger, fragmaxDir, sample, panddaDir, eventdf)

def export_models(logger, panddaDir, fragmaxDir, overwrite):
    sample_dict = get_sample_dict_for_export(logger, panddaDir)
    logger.info('starting sample export')
    for sample in sample_dict:
        logger.info('>>> {0!s} >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'.format(sample))
        ensemble = sample_dict[sample]
        pandda_model = os.path.join(panddaDir, 'processed_datasets', sample, 'modelled_structures', sample + '-pandda-model.pdb')
        if os.path.isfile(pandda_model):
            prepare_sample_refine_folder(logger, fragmaxDir, sample)
            if sample_refine_folder_is_empty(logger, fragmaxDir, sample):
                copy_files(logger, panddaDir, sample, fragmaxDir, ensemble)
            elif overwrite:
                backup_existing_files_folders_in_sample_dir(logger, fragmaxDir, sample)
                copy_files(logger, panddaDir, sample, fragmaxDir, ensemble)
            else:
                logger.warning("{0!s}: skipping sample; use -o option if you want to export...".format(sample))
        logger.info('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
    logger.info('finished exporting samples')


def main(argv):
    panddaDir = ''
    fragmaxDir = ''
    overwrite = False
    logger = processlib.init_logger('pandda_export.log')
    processlib.start_logging(logger, 'pandda_export.py')

    try:
        opts, args = getopt.getopt(argv,"p:f:ho",["panddadir=", "fragmaxdir=", "help", "overwrite"])
    except getopt.GetoptError:
        processlib.usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            processlib.usage()
            sys.exit(2)
        elif opt in ("-p", "--panddadir"):
            panddaDir = os.path.abspath(arg)
        elif opt in ("-f", "--fragmaxdir"):
            fragmaxDir = os.path.abspath(arg)
        elif opt in ("-o", "--overwrite"):
            overwrite = True

    if os.path.isdir(panddaDir) and os.path.isdir(fragmaxDir):
        export_models(logger, panddaDir, fragmaxDir, overwrite)
    else:
        logger.error('pandda directory and/or fragmax project folder do not exist: {0!s}'.format(panddaDir))

if __name__ == '__main__':
    main(sys.argv[1:])
