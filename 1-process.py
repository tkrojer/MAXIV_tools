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
            proc_list, match_found = processlib.retain_results_with_similar_ucvol_and_pg_as_ref_pdb(logger, proc_list, ref_dict)
            if not match_found:
                not_fitting_pipeline_list.append(sample)
                continue
#            proc_dict = processlib.retain_results_with_good_low_reso_rmerge(logger, proc_dict)
            proc_list = processlib.retain_results_with_good_low_reso_rmerge(logger, proc_list)
#            best, found_selected_pipeline = processlib.retain_results_which_fit_selection_criterion(logger, proc_dict, select_criterion)
            best, found_selected_pipeline = processlib.retain_results_which_fit_selection_criterion(logger, dal, proc_list, select_criterion)
            if best:
                processlib.link_process_results(logger, projectDir, sample, best)
                not_fitting_pipeline_list = processlib.check_if_best_result_is_from_select_pipeline(logger, sample, found_selected_pipeline, not_fitting_pipeline_list, select_criterion)
                processdb.unselected_autoprocessing_result(logger, dal, sample)
                processdb.set_selected_autoprocessing_result(logger, dal, sample, best)
            else:
                logger.error('None of MTZ files fulfilled the minimal requirements; check messages above')
#                processdb.unselected_autoprocessing_result(logger, dal, sample)
#                logger.warning('will select the one with the highest resolution')


        else:
            logger.error('could not find any MTZ or CIF in sample folder')
            processdb.select_last_dataset(logger, dal, sample)
            processlib.link_info_json_file(logger, projectDir, sample)
    processlib.report_not_fitting_pipelines(logger, not_fitting_pipeline_list, processDir, projectDir, fragmaxcsv)
    processlib.end_select_results(logger)


def parse_sample_folder(logger, sample_folder, projectDir, sample, proposal, session, pipelines,
                        status, protein, processDir, overwrite, beamline, dal, db_file, category, missing_dict, search_manual, n_manual):
    foundDataset = False
    foundMTZ = False
    for runs in sorted(glob.glob(os.path.join(sample_folder, '*'))):
        if not os.path.isdir(runs):
            continue
        if not os.listdir(runs):
            logger.error('process folder is empty: {0!s}'.format(runs))
            continue
#        run = runs.split('/')[9]
        run = runs.split('/')[len(runs.split('/'))-1]
        logger.info('checking run {0!s}'.format(run))

        if search_manual:
            logger.warning('looking for reprocessed fiiles only; assuming dataset table entry exists; skipping this part...')
        else:
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

#            logger.info('pipeline: {0!s}'.format(pipeline))
#            logger.info('projectDir: {0!s}'.format(projectDir))
#            logger.info('sample: {0!s}'.format(sample))
#            logger.info('run: {0!s}'.format(run))
#            logger.info('mtzpath: {0!s}'.format(mtzpath))

            if '_manual' in pipeline:
#                glob_string = os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run), mtzpath)
                glob_string = os.path.join(projectDir, '1-process', sample, '*', mtzpath)
                logger.info('glob_string: {0!s}'.format(glob_string))
            else:
                glob_string = os.path.join(sample_folder, run, mtzpath)

            for mtzfile in glob.glob(glob_string):
#                if '_manual' in pipeline:
                if search_manual:
                    logger.info('mtzfile manual: {0!s}'.format(mtzfile))
                    json_file = os.path.join(mtzfile[:mtzfile.rfind(pipeline)], "info.json")
                    if os.path.isfile(json_file):
                        logger.info("info.json exists in {0!s}".format(json_file))
                        proposal, session, run, collection_date = processlib.read_info_json_file(logger, json_file)
                        logger.info("copy json file to {0!s} folder".format(pipeline))
                        if not os.path.isfile(os.path.join(mtzfile[:mtzfile.rfind(pipeline)+len(pipeline)], "info.mtz")):
                            os.system("/bin/cp {0!s} {1!s}".format(json_file, mtzfile[:mtzfile.rfind(pipeline)+len(pipeline)]))
                    else:
                        logger.error("info.json does not exist in {0!s}; skipping...".format(json_file))
                        continue
                if processlib.process_files_for_run_pipeline_exist(logger, projectDir, sample, proposal, session, run,
                                                                   pipeline):
                    foundMTZ = True
#                    continue
                logger.info('found auto-processed MTZ file: ' + mtzfile)
                foundMTZ = True
                n_manual += 1
                status, logfile, ciffile, mtzfile, mrfana_ciffile = processlib.get_process_files(logger, mtzfile, projectDir, sample, proposal, session,
                                                      run, pipeline, collection_date,
                                                      mtz_extension, cif_extension, log_extension, status, mtz_unmerged)
                if not search_manual:
                    processlib.write_json_info_file(logger, projectDir, sample, collection_date, run, proposal, session,
                                                protein, status, master, pipeline)
                if ciffile or mrfana_ciffile:
                    d_xray_processing_table_dict = processdb.get_process_stats_from_mmcif_as_dict(logger, dal, ciffile, mtzfile,
                                                                                              logfile,
                                                                                              sample,
                                                                                              proposal, session, run, pipeline, projectDir, mrfana_ciffile)
                    if '_manual' in pipeline:
                        d_xray_processing_table_dict['automatic_processed'] = False
#                if os.path.isfile(db_file) and ciffile:
                if os.path.isfile(db_file) and d_xray_processing_table_dict:
                    if processdb.cif_exists(dal, ciffile) and not overwrite:
                        logger.info('cif file exisits in database and overwrite is False; skipping...')
                        continue

                    processdb.insert_into_xray_processing_table(logger, dal, d_xray_processing_table_dict, overwrite)
                    logger.info('finding highest resolution')
#                    print(d_xray_processing_table_dict)
#                    sys.exit()
                    processdb.assign_dataset_outcome(logger, dal, d_xray_processing_table_dict)

        if not search_manual:
            if not foundMTZ:
                if d_xray_dataset_table_dict:
                    processlib.write_json_info_file(logger, projectDir, sample, collection_date, run, proposal, session,
                                                    protein, status, master, None)

    if not foundDataset :
        if not search_manual:
            logger.warning('could not find any DATASET for sample, will create dummy entry in database...')
            processdb.create_dummy_dataset_entry(logger, dal, sample, proposal, session)
            missing_dict['dataset'].append([sample, proposal, session])
    if foundDataset and not foundMTZ:
        if not search_manual:
            missing_dict['mtz_file'].append([sample, proposal, session])
            logger.warning('could not find any MTZ file for sample!')
            status = processlib.get_status(logger, None, None, None, status)
            processlib.write_json_info_file(logger, projectDir, sample, collection_date, run, proposal, session,
                                            protein, status, master, '')
    logger.info('===================================================================================\n')
    return missing_dict, n_manual

def get_autoprocessing_results(logger, processDir, projectDir, fragmaxcsv, overwrite, dal, db_file, search_manual):
    sampleList = processlib.get_sample_list(logger, fragmaxcsv)
    missing_dict = {
        'dataset': [],
        'mtz_file': []
    }
    n_manual = 0
    if search_manual:
        folder = os.path.join(projectDir, "1-process")
    else:
        folder = processDir
        proposal, session, protein, beamline, category = processlib.get_proposal_and_session_and_protein(processDir)
        pipelines = processlib.get_processing_pipelines()

    for n, sample_folder in enumerate(sorted(glob.glob(os.path.join(folder, '*')))):
        sample = sample_folder.split('/')[len(sample_folder.split('/')) - 1]
        logger.info('current sample - {0!s}'.format(sample))
        if search_manual:
            logger.info("looking only for manually/ reprocessed data...")
            proposal, session, protein, beamline, category = processlib.get_none_proposal_and_session_and_protein()
            pipelines = processlib.get_manual_processing_pipelines()

        if sample in sampleList:
            logger.info('SUCCESS: found sample in summary csv file')
            processlib.create_sample_folder(logger, projectDir, sample)
            status = 'FAIL - no processing result'
            missing_dict, n_manual = parse_sample_folder(logger, sample_folder, projectDir, sample, proposal, session, pipelines,
                                status, protein, processDir, overwrite, beamline, dal, db_file, category, missing_dict, search_manual, n_manual)
        else:
            logger.warning('WARNING: cannot find sample in summary csv file')
            logger.info('===================================================================================\n')
    if search_manual:
        logger.info('found {0!s} MTZ files from manual processing'.format(n_manual))
    else:
        review_missing_datasets(logger, missing_dict, dal)
        processlib.reprocess_missing_datasets(logger, missing_dict, processDir, projectDir, fragmaxcsv)
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



def reprocess_datasets(logger, processDir, projectDir, reprocesscsv, overwrite, proc_dict, dal):
#    proposal, session, protein, beamline, category = processlib.get_proposal_and_session_and_protein(processDir)
    sampleList = processlib.get_sample_list(logger, reprocesscsv)
    n_jobs = 10  # hardcoded so that we don't hog the cluster
    now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    pipeline = proc_dict['pipeline'] + '_manual'
    script_dict = processlib.get_script_dict(pipeline, n_jobs, now)
    logger.info(script_dict)
    counter = 0
    for sample in sampleList:
        logger.info('current sample - {0!s}'.format(sample))
        master_files_runs = processdb.get_master_file_run_list(logger, dal, sample)
        processlib.create_sample_folder(logger, projectDir, sample)
        for item in master_files_runs:
            master_file = item[0]
            if not master_file:
                continue
            proposal, session, protein, beamline, category = processlib.get_proposal_and_session_and_protein(master_file)
            run = item[1]
            processlib.create_proposal_session_run_folder(logger, projectDir, sample, proposal, session, run)
            cont = processlib.create_pipeline_folder(logger, projectDir, sample, proposal, session, run, pipeline, overwrite)
            if cont:
                proc_folder = processlib.get_proc_folder(projectDir, sample, proposal, session, run, pipeline)
                script_dict = processlib.add_cmd_to_script_dict(logger, script_dict, counter, pipeline, proc_dict,
                                                    proc_folder, master_file, now)
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
    search_manual = False
    logger = processlib.init_logger('1-process.log')
    processlib.start_logging(logger, '1-process.py')

    try:
        opts, args = getopt.getopt(argv,"i:o:f:c:r:d:hsxm",["input=", "output=", "fragmax=", "criterion=", "reprocess=", "database=",
                                                       "help", "select", "overwrite", "manual"])
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
        elif opt in ("-m", "--manual"):
            search_manual = True

    processlib.report_parameters(logger, processDir, projectDir, fragmaxcsv, select, select_criterion, overwrite)
    checks_passed = processlib.run_checks(logger, processDir, projectDir, fragmaxcsv, select, select_criterion, db_file)

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
            reprocess_datasets(logger, processDir, projectDir, reprocesscsv, overwrite, proc_dict, dal)
        else:
            processlib.start_get_autoprocessing_results(logger)
            get_autoprocessing_results(logger, processDir, projectDir, fragmaxcsv, overwrite, dal, db_file, search_manual)
    else:
        logger.error('cannot continue; check error messages above and use -h option to get more information')
        logger.info('===================================================================================')

if __name__ == '__main__':
    main(sys.argv[1:])
