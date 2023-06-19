# Copyright (c) 2023, Tobias Krojer, MAX IV Laboratory
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

def get_datasets(project_directory, mtzin):
    mtz_list = []
    for mtz in sorted(glob.glob(os.path.join(project_directory, '*', mtzin))):
        mtz_list.append(mtz.split('/')[len(mtz.split('/'))-2])
    print('-> found {0!s} datasets'.format(len(mtz_list)))
    return mtz_list

def get_nproc():
    nproc = os.cpu_count()
    print('-> found {0!s} cpus in your system'.format(nproc))
    return nproc

def get_cmd_dict(nproc):
    cmd_dict = {}
    for i in range(nproc):
        cmd_dict['batch_{0!s}'.format(i)] = ''
    return cmd_dict

def make_cmd_dict(project_directory, mtzin, reference_pdb, reference_mtz, software):
    nproc = get_nproc()
    cmd_dict = get_cmd_dict(nproc)
    mtz_list = get_datasets(project_directory, mtzin)
    i = 0
    for sample in mtz_list:
        if n == nproc:
            n = 0
        cmd_dict['batch_{0!s}'.format(i)] += 'cd {0!s}\n'.format(os.path.join(project_directory,sample))
        cmd_dict['batch_{0!s}'.format(i)] += 'dimple {0!s} {1!s} {2!s} {3!s}\n'.format(mtzin,
                                                                                       reference_pdb,
                                                                                       reference_mtz,
                                                                                       software)
        n += 1
    return cmd_dict

def run_initial_refinement(project_directory, mtzin, reference_pdb, reference_mtz, software):
    cmd_dict = make_cmd_dict(project_directory, mtzin, reference_pdb, reference_mtz, software)
    for job in cmd_dict:
        print(cmd_dict[job])

def usage():
    print('hallo')

def main(argv):
    project_directory = ''
    overwrite = False
    software = 'dimple'
    reference_pdb = ''
    reference_mtz = ''
    mtzin = "process.mtz"

    try:
        opts, args = getopt.getopt(argv, "p:r:s:ho", ["project=", "reference=", "software=", "help", "overwrite"])
    except getopt.GetoptError:
        refinelib.usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit(2)
        elif opt in ("-p", "--project"):
            project_directory = os.path.abspath(arg)
        elif opt in ("-r", "--reference"):
            reference_pdb = os.path.abspath(arg)
        elif opt in ("-o", "--overwrite"):
            overwrite = True

    run_initial_refinement(project_directory, mtzin, reference_pdb, reference_mtz, software)

if __name__ == '__main__':
    main(sys.argv[1:])
