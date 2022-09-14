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

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'lib')))
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
        else:
            logger.error('could not find any MTZ or CIF in sample folder')
            processlib.link_info_json_file(logger, projectDir, sample)
    processlib.end_select_results(logger)


def parse_sample_folder(logger, sample_folder, projectDir, sample, proposal, session, pipelines,
                        status, protein, processDir, overwrite):
    foundMTZ = False
    for runs in glob.glob(os.path.join(sample_folder, '*')):
#        run = runs.split('/')[9]
        run = runs.split('/')[len(runs.split('/'))-1]
        logger.info('checking run {0!s}'.format(run))
        processlib.prepare_folders_and_files(logger, projectDir, sample, proposal, session, run, protein, processDir)
        collection_date, master = processlib.get_timestamp_from_master_file(sample_folder, run)
        for pipeline in pipelines:
            logger.info('checking {0!s} pipeline'.format(pipeline))
            mtzpath, mtz_extension, log_extension, cif_extension = processlib.get_pipeline_path(pipeline)
#            if processlib.process_files_for_run_pipeline_exist(logger, projectDir, sample, proposal, session, run, pipeline):
#                continue
#            logger.info('manual: ' + os.path.join(projectDir, '1-process', sample, '*', pipeline + '_*', mtz_extension))
            for mtzfile in glob.glob(os.path.join(sample_folder, '*', mtzpath)):
                if processlib.process_files_for_run_pipeline_exist(logger, projectDir, sample, proposal, session, run,
                                                                   pipeline):
                    foundMTZ = True
                    continue
                logger.info('found auto-processed MTZ file: ' + mtzfile)
                foundMTZ = True
                status = processlib.get_process_files(logger, mtzfile, projectDir, sample, proposal, session,
                                                      run, pipeline, collection_date,
                                                      mtz_extension, cif_extension, log_extension, status)
                processlib.write_json_info_file(logger, projectDir, sample, collection_date, run, proposal, session,
                                                protein, status, master, pipeline)
            # looking for manually processed datasets
            for mtzfile in glob.glob(os.path.join(projectDir, '1-process', sample, '*', pipeline + '_*', mtz_extension)):
 #               if processlib.process_files_for_run_pipeline_exist(logger, projectDir, sample, proposal, session, run,
 #                                                                  pipeline):
 #                   foundMTZ = True
 #                   continue
                logger.info('found manually processed MTZ file: ' + mtzfile)
                status = processlib.get_process_files(logger, mtzfile, projectDir, sample, proposal, session,
                                                      run, pipeline, collection_date,
                                                      mtz_extension, cif_extension, log_extension, status)
                processlib.write_json_info_file(logger, projectDir, sample, collection_date, run, proposal, session,
                                                protein, status, master, pipeline)


    if not foundMTZ:
        logger.warning('could not find any MTZ file for sample!')
        status = processlib.get_status(logger, None, None, None, status)
        processlib.write_json_info_file(logger, projectDir, sample, collection_date, run, proposal, session,
                                        protein, status, master, '')
    logger.info('===================================================================================\n')


def get_autoprocessing_results(logger, processDir, projectDir, fragmaxcsv, overwrite):
    sampleList = processlib.get_sample_list(logger, fragmaxcsv)
    proposal, session, protein = processlib.get_proposal_and_session_and_protein(processDir)
    pipelines = processlib.get_processing_pipelines()
    for n, sample_folder in enumerate(sorted(glob.glob(os.path.join(processDir, '*')))):
        sample = sample_folder.split('/')[len(sample_folder.split('/')) - 1]
        logger.info('current sample - {0!s}'.format(sample))
        if sample in sampleList:
            logger.info('SUCCESS: found sample in summary csv file')
            processlib.create_sample_folder(logger, projectDir, sample)
            status = 'FAIL - no processing result'
            parse_sample_folder(logger, sample_folder, projectDir, sample, proposal, session, pipelines,
                                status, protein, processDir, overwrite)
        else:
            logger.warning('WARNING: cannot find sample in summary csv file')
            logger.info('===================================================================================\n')
    processlib.end_get_autoprocessing_results(logger)


def reprocess_datasets(logger, processDir, projectDir, reprocesscsv, overwrite, proc_dict):
    proposal, session, protein = processlib.get_proposal_and_session_and_protein(processDir)
    sampleList = processlib.get_sample_list(logger, reprocesscsv)
    pipeline = proc_dict['pipeline'] + '_manual'
    n_jobs = 4  # hardcoded so that we don't hog the cluster
    script_dict = processlib.get_script_dict(pipeline, n_jobs)
    counter = 0
    for n, sample_folder in enumerate(sorted(glob.glob(os.path.join(processDir.replace('/process/', '/raw/'), '*')))):
        sample = sample_folder.split('/')[len(sample_folder.split('/')) - 1]
        if sample in sampleList:
            logger.info('current sample - {0!s}'.format(sample))
            processlib.create_sample_folder(logger, projectDir, sample)
            for master_file in glob.glob(os.path.join(sample_folder, '*_master.h5')):
                run = 'xds_' + master_file[master_file.rfind('/')+1:].replace('_master.h5', '') + '_1'
                processlib.create_proposal_session_run_folder(logger, projectDir, sample, proposal, session, run)
                processlib.create_pipeline_folder(logger, projectDir, sample, proposal, session, run, pipeline)
                proc_folder = processlib.get_proc_folder(projectDir, sample, proposal, session, run, pipeline)
                script_dict = processlib.add_cmd_to_script_dict(logger, script_dict, counter, pipeline, proc_dict,
                                                                proc_folder, master_file)
                counter += 1
                if counter == n_jobs:
                    counter = 0
    processlib.save_proc_scripts(logger, projectDir, script_dict)
    processlib.end_reprocessing(logger)


def main(argv):
    processDir = ''
    projectDir = ''
    fragmaxcsv = ''
    select = False
    select_criterion = 'resolution'
    reprocesscsv = ''
    overwrite = False
    logger = processlib.init_logger('1-process.log')
    processlib.start_logging(logger, '1-process.py')

    try:
        opts, args = getopt.getopt(argv,"i:o:f:c:r:hsx",["input=", "output=", "fragmax=", "crtierion=",
                                                       "help", "select", "overwrite", "reprocess="])
    except getopt.GetoptError:
        processlib.usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            processlib.usage()
            sys.exit(2)
        elif opt in ("-i", "--input"):
            processDir = os.path.abspath(arg)
        elif opt in ("-o", "--output"):
            projectDir = os.path.abspath(arg)
        elif opt in ("-f", "--fragmaxcsv"):
            fragmaxcsv = os.path.abspath(arg)
        elif opt in ("-s", "--select"):
            select = True
        elif opt in ("-c", "--crtierion"):
            select_criterion = arg
        elif opt in ("-x", "--overwrite"):
            overwrite = True
        elif opt in ("-r", "--reprocess"):
            reprocesscsv = os.path.abspath(arg)

    processlib.report_parameters(logger, processDir, projectDir, fragmaxcsv, select, select_criterion, overwrite)
    checks_passed = processlib.run_checks(logger, processDir, projectDir, fragmaxcsv, select, select_criterion)

    if checks_passed:
        processlib.check_if_to_continue(logger)
        if select:
            processlib.start_select_results(logger)
            select_results(logger, projectDir, select_criterion, overwrite)
        elif reprocesscsv:
            proc_dict = processlib.ask_for_spg_and_unit_cell(logger)
            processlib.start_reprocessing(logger)
            reprocess_datasets(logger, processDir, projectDir, reprocesscsv, overwrite, proc_dict)
        else:
            processlib.start_get_autoprocessing_results(logger)
            get_autoprocessing_results(logger, processDir, projectDir, fragmaxcsv, overwrite)
    else:
        logger.error('cannot continue; check error messages above and use -h option to get more information')
        logger.info('===================================================================================')

if __name__ == '__main__':
    main(sys.argv[1:])
