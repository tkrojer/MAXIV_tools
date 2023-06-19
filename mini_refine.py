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
import subprocess

def reference_fits_mtz(mtzfile, pdbDict):
    compatible = False
    try:
        import gemmi
        mtzDict = mtz_info(mtzfile)
        ucv_diff = abs((mtzDict['unitcell_volume'] - pdbDict['unitcell_volume'])) / pdbDict['unitcell_volume']
        if mtzDict['point_group'] == mtzDict['point_group'] and mtzDict['lattice'] == pdbDict['lattice'] and ucv_diff < 0.1:
            compatible = True
        else:
            print('ERROR: unit_cell volume, point group and lattice of reference and mtz file are too different')
    except ModuleNotFoundError:
        print('WARNING: cannot find gemmi module; cannot check if reference fits mtz file')
    return compatible

def mtz_info(mtzfile):
    mtzDict = {}
    mtz = gemmi.read_mtz_file(mtzfile)
    mtzDict['unitcell_volume'] = mtz.cell.volume
    mtzDict['point_group'] = mtz.spacegroup.point_group_hm()
    mtzDict['lattice'] = mtz.spacegroup.hm[0]
    return mtzDict

def pdb_info(pdbfile):
    pdbDict = {}
    try:
        import gemmi
        structure = gemmi.read_pdb(pdbfile)
        unitcell = structure.cell
        sym = gemmi.find_spacegroup_by_name(structure.spacegroup_hm)
        pdbDict['unitcell_volume'] = unitcell.volume
        mtzDict['point_group'] = sym.point_group_hm()
        pdbDict['lattice'] = sym.centring_type()
    except ModuleNotFoundError:
        print('WARNING: cannot find gemmi module; cannot check if reference fits mtz file')
    return pdbDict

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
#        cmd_dict['batch_{0!s}'.format(i)] = ''
        cmd_dict['batch_{0!s}'.format(i)] = 'bash -c "'
    return cmd_dict

def add_initial_refine_job(project_directory, mtzin, reference_pdb, reference_mtz, software, cmd_dict, i, sample):
#    script = 'bash -c "command1; command2; command3; " &'
    cmd_dict['batch_{0!s}'.format(i)] += 'cd {0!s};'.format(os.path.join(project_directory, sample))
    cmd_dict['batch_{0!s}'.format(i)] += 'dimple {0!s} {1!s} {2!s} {3!s};'.format(mtzin,
                                                                                   reference_pdb,
                                                                                   reference_mtz,
                                                                                   software)
    return cmd_dict

def make_cmd_dict(project_directory, mtzin, reference_pdb, reference_mtz, software):
    nproc = get_nproc()
    cmd_dict = get_cmd_dict(nproc)
    mtz_list = get_datasets(project_directory, mtzin)
    pdb_dict = pdb_info(reference_pdb)
    i = 0
    print('-> preparing initial refine script...')
    n_jobs = 0
    for sample in mtz_list:
        if i == nproc:
            i = 0
        if pdb_dict:
            if reference_mtz(os.path.join(project_directory,sample, mtzin), pdb_dict):
                cmd_dict = add_initial_refine_job(project_directory, mtzin, reference_pdb, reference_mtz, software,
                                                  cmd_dict, i, sample)
                i += 1
                n_jobs += 1
        else:
            cmd_dict = add_initial_refine_job(project_directory, mtzin, reference_pdb, reference_mtz, software,
                                              cmd_dict, i, sample)
            i += 1
            n_jobs += 1
    return cmd_dict, n_jobs

def submit(n_jobs):
    print('found {0!s} mtz files to process'.format(n_jobs))
    q = input("\n>>> Do you want to continue? (y/n) ")
    if not q.lower() == 'y':
        logger.info('you chose not to continue at this point; exciting program...')
        sys.exit(2)

def run_initial_refinement(project_directory, mtzin, reference_pdb, reference_mtz, software):
    cmd_dict, n_jobs = make_cmd_dict(project_directory, mtzin, reference_pdb, reference_mtz, software)
    submit(n_jobs)
    for job in cmd_dict:
        if cmd_dict[job]:
            script = cmd_dict[job] + ' " &'
            print(script)
#            subprocess.Popen([cmd_dict[job]])
#            os.spawnl(os.P_DETACH, cmd_dict[job])
#            os.system('{0!s} &'.format(cmd_dict[job]))

def usage():
    usage = (
        '\n'
        'usage:\n'
        'python mini_refine.py -d <data_dir> -p <reference_pdb>\n'
        '\n'
        'additional command line options:\n'
        '--data, -d\n'
        '    data directory\n'
        '--pdbref, -p\n'
        '    reference pdb file\n'
        '--mtzref, -r\n'
        '    reference mtz file\n'
        '--mtzin, -m\n'
        '    reference mtz file\n'
        '--overwrite, -o\n'
        '    flag to overwrite selected files\n'
        '\n'
        'Note: a) the script only works with python3\n'
        '      b) you need to install gemmi if you want to ensure that mtz and reference pdb file are similar'
    )
    print(usage)

def main(argv):
    project_directory = ''
    overwrite = False
    software = 'dimple'
    reference_pdb = ''
    reference_mtz = ''
    mtzin = "process.mtz"

    try:
        opts, args = getopt.getopt(argv, "d:p:r:m:ho", ["data=", "pdbref=", "mtzref=", "mtzin=",
                                                        "help", "overwrite"])
    except getopt.GetoptError:
        refinelib.usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit(2)
        elif opt in ("-d", "--data"):
            project_directory = os.path.abspath(arg)
        elif opt in ("-p", "--pdbref"):
            reference_pdb = os.path.abspath(arg)
        elif opt in ("-r", "--mtzref"):
            reference_mtz = os.path.abspath(arg)
        elif opt in ("-m", "--mtzin"):
            mtzin = os.path.abspath(arg)
        elif opt in ("-o", "--overwrite"):
            overwrite = True

    if not os.path.isdir(project_directory):
        print('ERROR: project directory does not exist; use -p flag to specify')
    elif not os.path.isfile(reference_pdb):
        print('ERROR: reference pdb file does not exist; use -r flag to specify')
    else:
        run_initial_refinement(project_directory, mtzin, reference_pdb, reference_mtz, software)

if __name__ == '__main__':
    main(sys.argv[1:])
