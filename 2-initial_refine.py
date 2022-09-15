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
from datetime import datetime

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'lib')))
import refinelib
import processlib


def run_initial_refinement(logger, projectDir, fragmaxcsv, software, overwrite):
    ref_dict = refinelib.get_reference_file_information(logger, projectDir)
    submitList = []
    counter = 0
    for l in open(fragmaxcsv):
        sample = l.split(',')[0]
        logger.info('current sample ' + sample)
        if refinelib.autoprocessing_files_exist(logger, projectDir, sample):
            mtzin = os.path.join(projectDir, '1-process', sample, 'process.mtz')
            mtzDict = processlib.mtz_info(mtzin)
            pdbref, mtzref = refinelib.suitable_reference_file_exists(logger, ref_dict, mtzDict)
            if pdbref:
                if refinelib.initial_refinement_exists(logger, projectDir, sample, software, overwrite):
                    now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
                    refinelib.prepare_script_for_init_refine(logger, projectDir, sample, mtzin, pdbref, mtzref,
                                                             now, submitList, counter, software)
                    counter += 1
    if submitList:
        logger.info('there are {0!s} {1!s} jobs to submit'.format(len(submitList), software))
        processlib.check_if_to_continue(logger)
        refinelib.submit_jobs_to_cluster(logger, projectDir, submitList)
    else:
        logger.warning('there are no jobs to submit; if this is unexpected, check messages above!')





def main(argv):
    projectDir = ''
    fragmaxcsv = ''
    overwrite = False
    linkrefine = False
    software = 'dimple'
    logger = processlib.init_logger('2-initial_refine.log')
    processlib.start_logging(logger, '2-initial_refine.py')
    try:
        opts, args = getopt.getopt(argv, "p:f:s:ho", ["project=", "fragmax=", "software=", "help", "overwrite"])
    except getopt.GetoptError:
        refinelib.usage()
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

#    processlib.report_parameters(logger, processDir, projectDir, fragmaxcsv, select, select_criterion, overwrite)
#    checks_passed = processlib.run_checks(logger, processDir, projectDir, fragmaxcsv, select, select_criterion)

    run_initial_refinement(logger, projectDir, fragmaxcsv, software, overwrite)

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
