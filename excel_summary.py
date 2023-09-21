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
import xlsxwriter
from PIL import Image
import getopt
import glob


sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'lib')))
import processlib
import summarylib


def parse_project_directory(logger, projectDir, fragmaxcsv, auxcsv):
    dataDict = summarylib.get_crystal_analysis_dict()
    statusDict = summarylib.get_status_dict()
    pgDict = summarylib.get_point_group_dict()
    pginfoDict = summarylib.get_point_group_info_dict()

    workbook = summarylib.init_workbook()
    n_samples = summarylib.get_n_samples(fragmaxcsv)
    n_autoproc_results = summarylib.n_autoproc_results(projectDir)
    summary_worksheet = summarylib.get_summary_worksheet(workbook, n_samples)
    details_worksheet = summarylib.get_details_worksheet(workbook, n_autoproc_results)
    pg_ucv_worksheet = summarylib.get_pg_ucv_worksheet(workbook)
    pg_rmergelow_worksheet = summarylib.get_pg_rmergelow_worksheet(workbook)

    chart_data = summarylib.get_chart_data_sheet(workbook)
    statusDict, df = summarylib.prepare_summary_worksheet(logger, workbook, summary_worksheet, projectDir, fragmaxcsv,
                                         dataDict, pgDict, pginfoDict, statusDict, auxcsv)
#    statusDict, df = summarylib.prepare_summary_worksheet_new(logger, workbook, summary_worksheet, projectDir, fragmaxcsv,
#                                         dataDict, pgDict, pginfoDict, statusDict, auxcsv)
    # this speeds up the process by a lot
#    summarylib.prepare_details_worksheet(logger, workbook, details_worksheet, projectDir, fragmaxcsv)
    summarylib.prepare_get_pg_ucv_worksheet(pg_ucv_worksheet)
    summarylib.prepare_get_pg_rmergelow_worksheet(pg_rmergelow_worksheet)
    summarylib.prepare_crystal_chart(workbook, chart_data, dataDict)
    summarylib.prepare_laue_group_chart(workbook, chart_data, pgDict)
    summarylib.prepare_status_chart(workbook, chart_data, statusDict)
    library_worksheet = summarylib.get_library_worksheet(workbook, df)
    workbook.close()
    df.to_csv('tsummary.csv', index=False)

def main(argv):
    projectDir = ''
    fragmaxcsv = ''
    auxcsv = ''
    logger = processlib.init_logger('excel_summary.log')
    processlib.start_logging(logger, 'excel_summary.py')

    try:
        opts, args = getopt.getopt(argv,"p:f:a:h",["project=", "fragmax=", "auxcsvn=", "help"])
    except getopt.GetoptError:
        processlib.usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            processlib.usage()
            sys.exit(2)
        elif opt in ("-p", "--project"):
            projectDir = os.path.abspath(arg)
        elif opt in ("-f", "--fragmaxcsv"):
            fragmaxcsv = os.path.abspath(arg)
        elif opt in ("-a", "--auxcsv"):
            auxcsv = os.path.abspath(arg)

#    processlib.report_parameters(logger, processDir, projectDir, fragmaxcsv, select, select_criterion, overwrite)
#    checks_passed = processlib.run_checks(logger, processDir, projectDir, fragmaxcsv, select, select_criterion)

    parse_project_directory(logger, projectDir, fragmaxcsv, auxcsv)

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
