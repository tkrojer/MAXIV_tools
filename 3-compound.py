# Copyright (c) 2022, Tobias Krojer, MAX IV Laboratory
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

import sys
import os
import getopt
import glob


sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'lib')))
import processlib
import summarylib
import compoundlib


def make_restraints(logger, projectDir, fragmaxcsv, software, overwrite):
    submitList = []
    for l in open(fragmaxcsv):
        sample = l.split(',')[0]
        if sample == 'SampleID' or sample == 'mounted_crystal_code':
            continue
        if not sample:
            print('no sample information in line; skipping...')
            continue
        logger.info('current sample ' + sample)
        cpdID = l.split(',')[4].replace(' ', '')
        smiles = l.split(',')[6].replace(' ', '')
        logger.info('CompoundID >{0!s}<, Smiles >{1!s}<'.format(cpdID, smiles))
        if smiles == '' or cpdID == '':
            logger.warning('CompoundID or Smiles field is empty; skipping...')
            continue
        compoundlib.create_sample_folder(logger, projectDir, sample)
        if cpdID:
            compoundlib.create_sample_folder(logger, projectDir, sample)
            if overwrite:
                compoundlib.empty_sample_folder(logger, projectDir, sample)
            compoundlib.make_compound_png(logger, projectDir, sample, cpdID, smiles)
            compoundlib.prepare_script_for_restrains_generation(logger, projectDir, sample, cpdID, smiles, software, submitList)
        else:
            logger.warning('sample does not have a compound ID/ smiles assigned; skipping')
    if submitList:
        logger.info('there are {0!s} {1!s} jobs to submit'.format(len(submitList), software))
        processlib.check_if_to_continue(logger)
        compoundlib.submit_jobs_to_cluster(logger, projectDir, submitList)
    else:
        logger.warning('there are no jobs to submit; if this is unexpected, check messages above!')


def make_links_to_pandda_folder(logger, projectDir, panddaDir):
    logger.info('making links from 3-compounds to {0!s}'.format(panddaDir))
    os.chdir(os.path.join(panddaDir, 'processed_datasets'))
    for sample_folder in glob.glob(os.path.join(panddaDir, 'processed_datasets', '*')):
        found_cif = False
        sample = sample_folder.split('/')[len(sample_folder.split('/'))-1]
        logger.info('current sample: ' + sample)
        for cif in glob.glob(os.path.join(projectDir, '3-compound', sample, '*.cif')):
            logger.info('found ligand cif file: ' + cif)
            found_cif = True
            break
        if found_cif:
            os.chdir(os.path.join(sample_folder, 'ligand_files'))
            os.system('/bin/rm *.cif')
            os.system('/bin/rm *.pdb')
#            os.system('ln -s {0!s}/*.pdb .'.format(os.path.join(projectDir, '3-compound', sample)))
#            os.system('ln -s {0!s}/*.cif .'.format(os.path.join(projectDir, '3-compound', sample)))
            os.system('/bin/cp {0!s}/*.pdb .'.format(os.path.join(projectDir, '3-compound', sample)))
            os.system('/bin/cp {0!s}/*.cif .'.format(os.path.join(projectDir, '3-compound', sample)))
        else:
            logger.warning('could not find ligand cif file')


def make_links_to_any_folder(logger, projectDir, anyDir):
    logger.info('making links from 3-compounds to {0!s}'.format(anyDir))
    os.chdir(anyDir)
    for sample_folder in glob.glob(os.path.join(anyDir, '*')):
        found_cif = False
        sample = sample_folder.split('/')[len(sample_folder.split('/'))-1]
        logger.info('current sample: ' + sample)
        for cif in glob.glob(os.path.join(projectDir, '3-compound', sample, '*.cif')):
            logger.info('found ligand cif file: ' + cif)
            found_cif = True
            break
        if found_cif:
            os.chdir(sample_folder)
            if not os.path.isdir('ligand_files'):
                os.mkdir('ligand_files')
            os.chdir('ligand_files')
            os.system('/bin/rm *.cif')
            os.system('/bin/rm *.pdb')
#            os.system('ln -s {0!s}/*.pdb .'.format(os.path.join(projectDir, '3-compound', sample)))
#            os.system('ln -s {0!s}/*.cif .'.format(os.path.join(projectDir, '3-compound', sample)))
            os.system('/bin/cp {0!s}/*.pdb .'.format(os.path.join(projectDir, '3-compound', sample)))
            os.system('/bin/cp {0!s}/*.cif .'.format(os.path.join(projectDir, '3-compound', sample)))
        else:
            logger.warning('could not find ligand cif file')

def run_rhofit(logger, projectDir):
    logger.info('running rhofit on dimple results')



def main(argv):
    projectDir = ''
    fragmaxcsv = ''
    anyDir = ''
    overwrite = False
    linkpandda = False
    software = 'acedrg'
    logger = processlib.init_logger('3-compound.log')
    processlib.start_logging(logger, '3-compound.py')


    try:
        opts, args = getopt.getopt(argv,"p:f:a:hol",["project=", "fragmax=", "anydir=", "help", "overwrite", "linkfiles"])
    except getopt.GetoptError:
        compoundlib.usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            compoundlib.usage()
            sys.exit(2)
        elif opt in ("-p", "--project"):
            projectDir = os.path.abspath(arg)
        elif opt in ("-f", "--fragmaxcsv"):
            fragmaxcsv = os.path.abspath(arg)
        elif opt in ("-o", "--overwrite"):
            overwrite = True
        elif opt in ("-l", "--linkfiles"):
            linkpandda = True
        elif opt in ("-a", "--anydir"):
            anyDir = os.path.abspath(arg)

#    processlib.report_parameters(logger, processDir, projectDir, fragmaxcsv, select, select_criterion, overwrite)
#    checks_passed = processlib.run_checks(logger, processDir, projectDir, fragmaxcsv, select, select_criterion)

    if linkpandda:
#        make_links_to_pandda_folder(logger, projectDir, panddaDir)
        make_links_to_any_folder(logger, projectDir, anyDir)
    else:
        make_restraints(logger, projectDir, fragmaxcsv, software, overwrite)

#    if checks_passed:
#        processlib.check_if_to_continue(logger)
#        if select:
#            processlib.start_select_results(logger)
#            select_results(logger, projectDir, select_criterion, overwrite)
#        else:
#            processlib.start_get_autoprocessing_results(logger)
#            get_autoprocessing_results(logger, processDir, projectDir, fragmaxcsv)
#    else:
#        logger.error('cannot continue; check error messages above and use -h option to get more information')
#        logger.info('===================================================================================')

if __name__ == '__main__':
    main(sys.argv[1:])
