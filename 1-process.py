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
from datetime import datetime

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'lib')))
import processlib
import processdb

sys.path.append('/data/staff/biomax/tobias/software/MAXIV_tools/lib')
from db import dal


def select_results(logger, projectDir, select_criterion, overwrite, processDir, fragmaxcsv, dal):
    logger.info('selecting auto-processing results based on {0!s}'.format(select_criterion))
    ref_dict = processlib.read_reference_pdb_files(logger, projectDir)
    not_fitting_pipeline_list = []
    for sample_folder in sorted(glob.glob(os.path.join(projectDir, '1-process', '*'))):
        proc_dict = {}
        os.chdir(sample_folder)
        sample = sample_folder.split('/')[len(sample_folder.split('/'))-1]
        logger.info('current sample - {0!s}'.format(sample))
        if processlib.skip_sample_if_already_selected(logger, projectDir, sample, sample_folder, overwrite):
            continue
#        for ciffile in glob.glob(os.path.join('*', '*', 'process.cif')):
#            proc_dict = processlib.read_data_collection_stats(logger, ciffile, proc_dict)
        proc_list = processdb.get_processing_results_for_sample(logger, dal, sample)
#        if proc_dict:
#        logger.warning('--> {0!s}'.format(proc_list))
        if proc_list:
#            proc_dict = processlib.retain_results_with_similar_ucvol_and_pg_as_ref_pdb(logger, proc_dict, ref_dict)
            proc_list = processlib.retain_results_with_similar_ucvol_and_pg_as_ref_pdb(logger, proc_list, ref_dict)
#            proc_dict = processlib.retain_results_with_good_low_reso_rmerge(logger, proc_dict)
            proc_list = processlib.retain_results_with_good_low_reso_rmerge(logger, proc_list)
#            best, found_selected_pipeline = processlib.retain_results_which_fit_selection_criterion(logger, proc_dict, select_criterion)
            best, found_selected_pipeline = processlib.retain_results_which_fit_selection_criterion(logger, proc_list, select_criterion)
            if best:
                processlib.link_process_results(logger, projectDir, sample, best)
                not_fitting_pipeline_list = processlib.check_if_best_result_is_from_select_pipeline(logger, sample, found_selected_pipeline, not_fitting_pipeline_list, select_criterion)
                processdb.unselected_autoprocessing_result(logger, dal, sample)
                processdb.set_selected_autoprocessing_result(logger, dal, sample, best)
            else:
                logger.error('None of MTZ files fulfilled the minimal requirements; check messages aboove')
        else:
            logger.error('could not find any MTZ or CIF in sample folder')
            processlib.link_info_json_file(logger, projectDir, sample)
    processlib.report_not_fitting_pipelines(logger, not_fitting_pipeline_list, processDir, projectDir, fragmaxcsv)
    processlib.end_select_results(logger)


def parse_sample_folder(logger, sample_folder, projectDir, sample, proposal, session, pipelines,
                        status, protein, processDir, overwrite, beamline, dal, db_file, category, missing_dict):
    foundDataset = False
    foundMTZ = False
    for runs in sorted(glob.glob(os.path.join(sample_folder, '*'))):
#        run = runs.split('/')[9]
        run = runs.split('/')[len(runs.split('/'))-1]
        logger.info('checking run {0!s}'.format(run))



        dozor_plot, crystal_snapshot_list = processlib.prepare_folders_and_files(logger, projectDir, sample, proposal, session, run, protein, processDir, category, beamline)
        collection_date, master, create_date = processlib.get_timestamp_from_master_file(logger, sample_folder, run)

        d_xray_dataset_table_dict, foundDataset = processdb.get_d_xray_dataset_table_dict(logger, dal, sample, proposal, session, beamline,
                                                                            run, create_date, master, dozor_plot, crystal_snapshot_list, foundDataset)

        if os.path.isfile(db_file):
            processdb.insert_into_xray_dataset_table(logger, dal, d_xray_dataset_table_dict)

        if not master:  # this may happen if there is a run folder, but without image files
            continue
        for pipeline in pipelines:
            logger.info('checking {0!s} pipeline'.format(pipeline))
            mtzpath, mtz_extension, log_extension, cif_extension, mtz_unmerged = processlib.get_pipeline_path(pipeline)
            for mtzfile in glob.glob(os.path.join(sample_folder, run, mtzpath)):
                if processlib.process_files_for_run_pipeline_exist(logger, projectDir, sample, proposal, session, run,
                                                                   pipeline):
                    foundMTZ = True
                    continue
                logger.info('found auto-processed MTZ file: ' + mtzfile)
                foundMTZ = True
                status, logfile, ciffile, mtzfile = processlib.get_process_files(logger, mtzfile, projectDir, sample, proposal, session,
                                                      run, pipeline, collection_date,
                                                      mtz_extension, cif_extension, log_extension, status, mtz_unmerged)
                processlib.write_json_info_file(logger, projectDir, sample, collection_date, run, proposal, session,
                                                protein, status, master, pipeline)
                if ciffile:
                    d_xray_processing_table_dict = processdb.get_process_stats_from_mmcif_as_dict(logger, dal, ciffile, mtzfile,
                                                                                              logfile,
                                                                                              sample,
                                                                                              proposal, session, run)
                if os.path.isfile(db_file) and ciffile:
                    processdb.insert_into_xray_processing_table(logger, dal, d_xray_processing_table_dict)
                    logger.info('finding highest resolution')
                    processdb.assign_dataset_outcome(logger, dal, sample, d_xray_processing_table_dict)

    #            # looking for manually processed datasets
#            manual = pipeline
#            if pipeline == 'staraniso':
#                manual = 'autoproc'
#            for mtzfile in glob.glob(os.path.join(projectDir, '1-process', sample, '*', manual + '_*', mtz_extension)):
#                manual_pipeline = processlib.get_manual_pipeline_name(logger, manual, mtzfile)
#                if pipeline == 'staraniso':
#                    manual_pipeline = manual_pipeline.replace('autoproc', 'staraniso')
#
#                if processlib.process_files_for_run_pipeline_exist(logger, projectDir, sample, proposal, session, run,
#                                                                   manual_pipeline):
#                    foundMTZ = True
#                    continue
#                logger.info('found manually processed MTZ file: ' + mtzfile)
#                status, logfile = processlib.get_process_files(logger, mtzfile, projectDir, sample, proposal, session,
#                                                      run, manual_pipeline, collection_date,
#                                                      mtz_extension, cif_extension, log_extension, status)
#                processlib.write_json_info_file(logger, projectDir, sample, collection_date, run, proposal, session,
#                                                protein, status, master, manual_pipeline)

#    if foundDataset and foundMTZ:
#        logger.info('finding highest resolution')
#        processdb.assign_dataset_outcome(logger, dal, sample)


    if not foundDataset :
        logger.warning('could not find any DATASET for sample, will create dummy entry in database...')
        processdb.create_dummy_dataset_entry(logger, dal, sample, proposal, session)
        missing_dict['dataset'].append([sample, proposal, session])
    if foundDataset and not foundMTZ:
        missing_dict['mtz_file'].append([sample, proposal, session])
        logger.warning('could not find any MTZ file for sample!')
        status = processlib.get_status(logger, None, None, None, status)
        processlib.write_json_info_file(logger, projectDir, sample, collection_date, run, proposal, session,
                                        protein, status, master, '')
    logger.info('===================================================================================\n')
    return missing_dict

def get_autoprocessing_results(logger, processDir, projectDir, fragmaxcsv, overwrite, dal, db_file):
    sampleList = processlib.get_sample_list(logger, fragmaxcsv)
    proposal, session, protein, beamline, category = processlib.get_proposal_and_session_and_protein(processDir)
    pipelines = processlib.get_processing_pipelines()
    missing_dict = {
        'dataset': [],
        'mtz_file': []
    }
    for n, sample_folder in enumerate(sorted(glob.glob(os.path.join(processDir, '*')))):
        sample = sample_folder.split('/')[len(sample_folder.split('/')) - 1]
        logger.info('current sample - {0!s}'.format(sample))
        if sample in sampleList:
            logger.info('SUCCESS: found sample in summary csv file')
            processlib.create_sample_folder(logger, projectDir, sample)
            status = 'FAIL - no processing result'
            missing_dict = parse_sample_folder(logger, sample_folder, projectDir, sample, proposal, session, pipelines,
                                status, protein, processDir, overwrite, beamline, dal, db_file, category, missing_dict)
        else:
            logger.warning('WARNING: cannot find sample in summary csv file')
            logger.info('===================================================================================\n')
    review_missing_datasets(logger, missing_dict, dal)
    processlib.end_get_autoprocessing_results(logger)


def review_missing_datasets(logger, missing_dict, dal):
    logger.warning('the following samples have either no dataset collected or the dataset did not result in a MTZ file:')
    logger.info('missing datasets:')
    for i in missing_dict['dataset']:
        logger.info(' --> {0!s}'.format(i[0]))
    logger.info('missing MTZ file:')
    for i in missing_dict['mtz_file']:
        logger.info(' --> {0!s}'.format(i[0]))
    processlib.check_if_to_annotate_and_reprocess(logger, missing_dict, dal)



def reprocess_datasets(logger, processDir, projectDir, reprocesscsv, overwrite, proc_dict):
    proposal, session, protein, beamline, category = processlib.get_proposal_and_session_and_protein(processDir)
    sampleList = processlib.get_sample_list(logger, reprocesscsv)
    n_jobs = 4  # hardcoded so that we don't hog the cluster
    now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    pipeline = proc_dict['pipeline'] + '_' + now
    script_dict = processlib.get_script_dict(pipeline, n_jobs)
    counter = 0
    for sample in sampleList:
        logger.info('current sample - {0!s}'.format(sample))
        master_files_runs = processdb.get_master_file_run_list(logger, dal, sample, proposal, session)
        for item in master_files_runs:
            master_file = item[0]
            run = item[1]



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
    db_file = ""
    logger = processlib.init_logger('1-process.log')
    processlib.start_logging(logger, '1-process.py')

    try:
        opts, args = getopt.getopt(argv,"i:o:f:c:r:d:hsx",["input=", "output=", "fragmax=", "crtierion=", "database=",
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
        elif opt in ("-c", "--criterion"):
            select_criterion = arg
        elif opt in ("-x", "--overwrite"):
            overwrite = True
        elif opt in ("-r", "--reprocess"):
            reprocesscsv = os.path.abspath(arg)
        elif opt in ("-d", "--database"):
            db_file = os.path.abspath(arg)

    processlib.report_parameters(logger, processDir, projectDir, fragmaxcsv, select, select_criterion, overwrite)
    checks_passed = processlib.run_checks(logger, processDir, projectDir, fragmaxcsv, select, select_criterion)

    if checks_passed:
        processlib.check_if_to_continue(logger)
        logger.info('initializing database: {0!s}'.format(db_file))
        dal.db_init(db_file)
        if select:
            processlib.start_select_results(logger)
            select_results(logger, projectDir, select_criterion, overwrite, processDir, fragmaxcsv, dal)
        elif reprocesscsv:
            proc_dict = processlib.ask_for_spg_and_unit_cell(logger)
            processlib.start_reprocessing(logger)
            reprocess_datasets(logger, processDir, projectDir, reprocesscsv, overwrite, proc_dict)
        else:
            processlib.start_get_autoprocessing_results(logger)
            get_autoprocessing_results(logger, processDir, projectDir, fragmaxcsv, overwrite, dal, db_file)
    else:
        logger.error('cannot continue; check error messages above and use -h option to get more information')
        logger.info('===================================================================================')

if __name__ == '__main__':
    main(sys.argv[1:])
