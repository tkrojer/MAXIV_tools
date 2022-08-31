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


import getopt
import glob
import sys
import os
import gemmi
import time
import json
import logging
import math


sys.path.append('/data/staff/biomax/tobias/script/lib')
import processlib


def select_results(logger, projectDir, select_criterion, overwrite):
    logger.info('selecting auto-processing results based on {0!s}'.format(select_criterion))
    ref_dict = processlib.read_reference_pdb_files(logger, projectDir)
    for sample_folder in glob.glob(os.path.join(projectDir, '1-process', '*')):
        proc_dict = {}
        os.chdir(sample_folder)
        sample = sample_folder.split('/')[len(sample_folder.split('/'))-1]
        logger.info('current sample - {0!s}'.format(sample))
        if processlib.skip_sample_if_already_selected(logger, projectDir, sample, sample_folder, overwrite):
            continue
        for ciffile in glob.glob(os.path.join('*', '*', 'process.cif')):
            proc_dict = processlib.read_data_collection_stats(logger, ciffile, proc_dict)
        if proc_dict:
            proc_dict = processlib.retain_results_with_similar_ucvol_and_pg_as_ref_pdb(logger, proc_dict, ref_dict)
            proc_dict = processlib.retain_results_with_good_low_reso_rmerge(logger, proc_dict)
            best = processlib.retain_results_which_fit_selection_criterion(logger, proc_dict, select_criterion)
            processlib.link_process_results(logger, projectDir, sample, best)
        logger.error('could not find any MTZ or CIF in sample folder')
    processlib.end_select_results(logger)








def get_autoprocessing_results(logger, processDir, projectDir, fragmaxcsv):
    sampleList = processlib.get_sample_list(fragmaxcsv)
    proposal, session, protein = processlib.get_proposal_and_session_and_protein(processDir)
    pipelines = processlib.get_processing_pipelines()
    for n, sample_folder in enumerate(sorted(glob.glob(os.path.join(processDir, '*')))):
        sample = sample_folder.split('/')[len(sample_folder.split('/')) - 1]
        logger.info('current sample - {0!s}'.format(sample))
        if sample in sampleList:
            logger.info('SUCCESS: found sample in summary csv file')
            processlib.create_sample_folder(projectDir, sample)
            foundMTZ = False
            status = 'FAIL - no processing result'
            for runs in glob.glob(os.path.join(sample_folder, '*')):
                run = runs.split('/')[9]
                logger.info('checking run {0!s}'.format(run))
                processlib.prepare_folders_and_files(logger, projectDir, sample, proposal, session, run)
                collection_date, master = processlib.get_timestamp_from_master_file(sample_folder, run)
                for pipeline in pipelines:
                    logger.info('checking {0!s} pipeline'.format(pipeline))
                    mtzpath, mtz_extension, log_extension, cif_extension = processlib.get_pipeline_path(pipeline)
                    for mtzfile in glob.glob(os.path.join(sample_folder, '*', mtzpath)):
                        logger.info('found MTZ file')
                        foundMTZ = True
                        status = processlib.get_process_files(logger, mtzfile, projectDir, sample, proposal, session,
                                                              run, pipeline, collection_date,
                                                              mtz_extension, cif_extension, log_extension)
                if not foundMTZ:
                    logger.warning('could not find any MTZ file!')
                    status = processlib.get_status(None, status)
                processlib.write_json_info_file(logger, projectDir, sample, collection_date, run, proposal, session,
                                                protein, status, master)
            logger.info('\n===================================================================================\n')

        else:
            logger.warning('WARNING: cannot find sample in summary csv file')
    processlib.end_get_autoprocessing_results(loger)


def check_if_process_directory_exists(processDir):
    passed = True
    print('-> checking process directory: {0!s}'.format(processDir))
    if os.path.isdir(processDir):
        print('OK: process directory exists')
    else:
        print('ERROR: process directory does not exisit')
        passed = False
    return passed

def check_if_project_directory_exists(projectDir, passed):
    print('-> checking project directory: {0!s}'.format(projectDir))
    if os.path.isdir(projectDir):
        print('OK: project directory exists')
    else:
        print('ERROR: project directory does not exisit')
        passed = False
    return passed


#def check_if_pdbdir_directory_exists(pdbDir, passed):
#    print('-> checking PDB directory: {0!s}'.format(pdbDir))
#    if os.path.isdir(pdbDir):
#        print('OK: PDB directory exists')
#    else:
#        print('ERROR: PDB directory does not exisit')
#        passed = False
#    return passed

#def check_pdbfiles_in_pdbdir(pdbDir, passed):
#    print('-> looking for PDB files in PDB directory...')
#    found = False
#    for pdb in glob.glob(os.path.join(pdbDir, '*.pdb')):
#        pdbFile = pdb.split('/')[len(pdb.split('/'))-1]
#        print('OK: found {0!s}'.format(pdbFile))
#        foundCryst = False
#        for line in open(pdb):
#            if line.startswith('CRYST'):
#                print('OK: found CRYST card: {0!s}'.format(line.replace('\n', '')))
#                foundCryst = True
#                found = True
#        if not foundCryst:
#            print('WARNING: {0!s} does not seem to contain a CRYST card'.format(pdbFile))
#    if found:
#        print('OK: found at least one PDB file with a CRYST card')
#    else:
#        print('ERROR: did not find a PDB file with a valid CRYST card')
#        passed = False
#    return passed


def check_process_pipeline_option(process_pipeline, passed):
    print('-> checking process pipeline option: {0!s}'.format(process_pipeline))
    supported_options = ['dials', 'autoproc', 'staraniso', 'fastdp']
    if process_pipeline in supported_options:
        print('OK: option exists')
    else:
        print('ERROR: option does not exist')
        passed = False
    return passed


def check_refine_pipeline_option(refine_pipeline, passed):
    print('-> checking refine pipeline option: {0!s}'.format(refine_pipeline))
    supported_options = ['dimple', 'pipedream']
    if refine_pipeline in supported_options:
        print('OK: option exists')
    else:
        print('ERROR: option does not exist')
        passed = False
    return passed


def run_checks(processDir, projectDir, process_pipeline, refine_pipeline):
    print('>>> checking input file and command line options')
    passed = check_if_process_directory_exists(processDir)
    passed = check_if_project_directory_exists(projectDir, passed)
#    passed = check_if_pdbdir_directory_exists(pdbDir, passed)
#    passed = check_pdbfiles_in_pdbdir(pdbDir, passed)
    passed = check_process_pipeline_option(process_pipeline, passed)
    passed = check_refine_pipeline_option(refine_pipeline, passed)
    return passed


def usage():
    usage = (
        '\n'
        'usage:\n'
        'ccp4-python run_initial_refinement.py -i <process_dir> -o <project_dir>\n'
        '\n'
        'additional command line options:\n'
        '--input, -i\n'
        '    process directory\n'
        '--output, -o\n'
        '    project directory\n'
        '--autoproc, -a\n'
        '    auto-processing pipeline (e.g. autooproc)\n'
        '--refine, -r\n'
        '    initial refinement pipeline (e.g. dimple)\n'
        '--overwrite, -o\n'
        '    flag to overwrite files\n'
        '--analyse, -y\n'
        '    flag to analyse process directory, without file operations\n'
        '--create, -c\n'
        '    create subfolders in project directory (Note: run this first)\n'
    )
    print(usage)


def check_if_to_continue():
    if sys.version[0] == '2':
        q = raw_input("\n>>> Do you want to continue? (y/n) ")
    else:
        q = input("\n>>> Do you want to continue? (y/n)")
    if not q.lower() == 'y':
        print('\n>>> exciting program...')
        sys.exit(2)










def main(argv):
    cwd = os.getcwd()
    processDir = None
    projectDir = None
    fragmaxcsv = None
    select = False
    select_criterion = 'resolution'
    logger = processlib.init_logger()

    processDir.start_info(logger)


    show_title()

    try:
        opts, args = getopt.getopt(argv,"i:o:f:hs",["input=", "output=", "fragmax=", "help", "select"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-i", "--input"):
            processDir = os.path.abspath(arg)
        elif opt in ("-o", "--output"):
            projectDir = os.path.abspath(arg)
        elif opt in ("-f", "--fragmax"):
            fragmaxcsv = os.path.abspath(arg)
        elif opt in ("-s", "--select"):
            select = True
        elif opt in ("-c", "--crtierion"):
            select_criterion = arg

#    if createSubdirectories:
#        print("WARNING: will ignore all other options and only create subfolders in project directory")
#        prepare_directory_structure(projectDir)
#        sys.exit("done")

#    checks_passed = run_checks(processDir, projectDir, process_pipeline, refine_pipeline)

#    if checks_passed:
#        check_if_to_continue()
##        prepare_directory_structure(projectDir)
    if select:
        processlib.start_select_results(logger)
        select_results(logger, projectDir, select_criterion)
    else:
        processlib.start_get_autoprocessing_results(logger)
        get_autoprocessing_results(logger, processDir, projectDir, fragmaxcsv)
#    else:
#        print('something is wrong, please check comments above and use -h option for more information')


if __name__ == '__main__':
    main(sys.argv[1:])
