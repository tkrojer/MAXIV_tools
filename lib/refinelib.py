import glob
import os
import sys

import gemmi


def usage():
    usage = (
        '\n'
        'usage:\n'
        'python 2-initial_refine.py -p <project_dir> -f <fragmax_csv_file>\n'
        '\n'
        'additional command line options:\n'
        '--project, -p\n'
        '    project directory\n'
        '--fragmax, -f\n'
        '    fragmax summary csv file\n'
        '--overwrite, -o\n'
        '    flag to overwrite selected files\n'
        '--software, -s\n'
        '    initial refinement pipeline (dimple [default], pipedreamo)\n'
    )
    print(usage)


def create_sample_folder(logger, projectDir, sample):
    os.chdir(os.path.join(projectDir, '2-initial_refine'))
    if os.path.isdir(sample):
        logger.warning('sample folder exists')
    else:
        logger.info('creating folder {0!s} in 2-initial_refine')
        os.mkdir(sample)


def get_reference_file_information(logger, projectDir):
    logger.info('looking for reference files in {0!s}/0-model...'.format(projectDir))
    ref_dict = {}
    ciffile = None
    for pdbfile in glob.glob(os.path.join(projectDir, '0-model', '*pdb')):
        if os.path.isfile(pdbfile.replace('.pdb', '.cif')):
            ciffile = pdbfile.replace('.pdb', '.cif')
        pdb = pdbfile[pdbfile.rfind('/')+1:]
        mtz = ''
        structure = gemmi.read_pdb(pdbfile)
        unitcell = structure.cell
        unitcell_volume = unitcell.volume
        sym = gemmi.find_spacegroup_by_name(structure.spacegroup_hm)
        point_group = sym.point_group_hm()
        lattice = sym.centring_type()
        if os.path.isfile(pdbfile.replace('.pdb', '.mtz')):
            mtz = pdb.replace('.pdb', '.mtz')
        logger.info('found {0!s} - pg {1!s} - lattice {2!s} - uc_vol {3!s} - ref_mtz {4!s}'.format(
            pdb, point_group, lattice, unitcell_volume, mtz
        ))
        ref_dict[pdbfile] = [point_group, lattice, unitcell_volume, mtz, ciffile]
    if not ref_dict:
        logger.error('could not find any PDB file in folder')
    return ref_dict


def autoprocessing_files_exist(logger, projectDir, sample):
    processmtz_exists = False
    logger.info('checking if auto-processing files exisit...')
    if os.path.isfile(os.path.join(projectDir, '1-process', sample, 'process.mtz')):
        logger.info('process.mtz file exists in "1-process/{0!s}" folder'.format(sample))
        processmtz_exists = True
    else:
        logger.error('process.mtz file does not exist in "1-process/{0!s}" folder; skipping sample...'.format(sample))
    return processmtz_exists


def suitable_reference_file_exists(logger, ref_dict, mtzDict):
    logger.info('looking for suitable input PDB file...')
    pdbref = None
    cifref = None
    mtzref = ''
    pgr_mtz = mtzDict['point_group']
    ucv_mtz = float(mtzDict['unitcell_volume'])
    lat_mtz = mtzDict['lattice']
    for pdb in ref_dict:
        pgr_pdb = ref_dict[pdb][0]
        ucv_pdb = float(ref_dict[pdb][2])
        lat_pdb = ref_dict[pdb][1]
        mtzref = ref_dict[pdb][3]
        if ref_dict[pdb][4]:
            cifref = ref_dict[pdb][4]
        ucv_diff = abs((ucv_mtz - ucv_pdb)) / ucv_pdb
        if pgr_mtz == pgr_pdb and lat_mtz == lat_pdb and ucv_diff < 0.1:
            logger.info('lattice, point group and unit cell volume of MTZ file matches {0!s}'.format(pdb))
            pdbref = pdb
            break
    if not pdbref:
        logger.error('could not find a suitable reference PDB file for MTZ file')
    return pdbref, mtzref, cifref


def initial_refinement_exists(logger, projectDir, sample, software, overwrite):
    logger.info('checking if initial refinement folder for {0!s} exists...'.format(software))
    continue_refinement = True
    if os.path.isdir(os.path.join(projectDir, '2-initial_refine', sample, software)):
        logger.warning('{0!s} folder exists'.format(software))
        if overwrite:
            logger.warning('"overwrite" option selected; will remove folder...')
            os.chdir(os.path.join(projectDir, '2-initial_refine', sample))
#            os.system('/bin/rm -fr ' + software)
        else:
            logger.error('you need to select "overwrite" option if you want to run init refinement; skipping sample...')
            continue_refinement = False
    return continue_refinement


def prepare_script_for_init_refine(logger, projectDir, sample, mtzin, pdbref, mtzref, now, submitList, counter, software, cifref, cpd_id):
    logger.info('preparing script for initial refinement')
    cmd = maxiv_header(software)
    cmd += modules_to_load(software)
    cmd += init_refine_cmd(logger, software, projectDir, sample, mtzin, pdbref, mtzref, cifref, cpd_id)
    os.chdir(os.path.join(projectDir, 'tmp'))
    submitList.append('{0!s}_{1!s}_{2!s}.sh'.format(software, now, counter))
    f = open('{0!s}_{1!s}_{2!s}.sh'.format(software, now, counter), 'w')
    f.write(cmd)
    f.close()


def maxiv_header(software):
    header = (
        '#!/bin/bash\n'
        '#SBATCH --time=10:00:00\n'
        '#SBATCH --job-name={0!s}\n'.format(software) +
        '#SBATCH --cpus-per-task=1\n'
    )
    return header


def modules_to_load(software):
    if software == 'dimple':
        module = 'module load gopresto CCP4\n'
    elif software == 'pipedream':
        module = 'module load gopresto BUSTER\n'
    elif software == 'phenix':
        module = 'module load gopresto Phenix\n'
    return module


def init_refine_cmd(logger, software, projectDir, sample, mtzin, pdbref, mtzref, cifref, cpd_id):
    cmd = 'cd {0!s}\n'.format(os.path.join(projectDir, '2-initial_refine', sample))
    if software == 'dimple':
        cmd += 'dimple {0!s} {1!s} {2!s} {3!s}\n'.format(mtzin, pdbref, mtzref, software)
    elif software == 'pipedream':
        if cifref:
            lig = "-l {0!s}".format(cifref)
        else:
            lig = ""
        rhofit = ""
        logger.info('looking for {0!s}'.format(os.path.join(projectDir, '3-compound', sample, cpd_id + '.cif')))
        if os.path.isdir(os.path.join(projectDir, '3-compound', sample, cpd_id + '.cif')):
            rhofit = "-rhofit {0!s}".format(os.path.join(projectDir, '3-compound', sample, cpd_id + '.cif'))
            cmd += 'pipedream -xyzin {0!s} -hklin {1!s} -nofreeref -nolmr -d {2!s} {3!s} {4!s}'.format(pdbref, mtzin,
                                                                                                       software, lig,
                                                                                                       rhofit)
        else:
            rhofit = ""
            cmd += ""
    elif software == 'phenix':
        cmd += ''
    return cmd


def submit_jobs_to_cluster(logger, projectDir, submitList):
    logger.info('submitting jobs to MAXIV cluster...')
    os.chdir(os.path.join(projectDir, 'tmp'))
    for script in submitList:
        logger.info('submitting ' + script)
#        os.system('sbatch ' + script)


#def structure_cif_info(cif):
#    cifDict = {
#        'ls_d_res_high': '',
#        'ls_R_factor_R_work': '',
#        'ls_R_factor_R_free': '',
#        'initref_software': '',
#        'initref_spacegroup': '',
#        'r_bond_refined_d': '',
#        'r_angle_refined_deg': ''
#    }
#    doc = gemmi.cif.read_file(cif)
#    for block in doc:
#        if block.find_pair('_refine.ls_d_res_high'):
#            cifDict['ls_d_res_high'] = str(block.find_pair('_refine.ls_d_res_high')[1])
#        if block.find_pair('_refine.ls_R_factor_R_work'):
#            cifDict['ls_R_factor_R_work'] = str(block.find_pair('_refine.ls_R_factor_R_work')[1])
#        if block.find_pair('_refine.ls_R_factor_R_free'):
#            cifDict['ls_R_factor_R_free'] = str(block.find_pair('_refine.ls_R_factor_R_free')[1])
#        if block.find_pair('_software.name'):
#            cifDict['initref_software'] = str(block.find_pair('_software.name')[1])
#        if block.find_pair('_symmetry.space_group_name_H-M'):
#            cifDict['initref_spacegroup'] = str(block.find_pair('_symmetry.space_group_name_H-M')[1])
#        if block.find_loop('_refine_ls_restr.type'):
#            table = block.find('_refine_ls_restr.', ['type', 'dev_ideal'])
#            cifDict['r_bond_refined_d'] = str(list(table.find_row('r_bond_refined_d'))[1])
#            cifDict['r_angle_refined_deg'] = str(list(table.find_row('r_angle_refined_deg'))[1])
#    return cifDict


def find_blobs(mtzfile, pdbfile):
    blobList = []
    mtz = gemmi.read_mtz_file(mtzfile)
    grid = mtz.transform_f_phi_to_map('FWT', 'PHWT', sample_rate=3)
    st = gemmi.read_structure(pdbfile)
    st.remove_waters()
    masker = gemmi.SolventMasker(gemmi.AtomicRadiiSet.Constant, 1.75)
    masker.set_to_zero(grid, st[0])
    blobs = gemmi.find_blobs_by_flood_fill(grid, cutoff=0.6, min_volume=8, min_score=0, min_peak=0)
    for blob in blobs:
        blobList.append(blob.volume)
    return blobList

"""
            initpdb, initmtz, initcif, freecif = processlib.get_refinement_files(logger, projectDir, sample, software)
            os.chdir(os.path.join(projectDir, '2-initial_refine', sample))
            if os.path.isfile(os.path.join(software, 'final.pdb')):
                os.system('ln -s {0!s} init.pdb'.format(os.path.join(software, 'final.pdb')))
            if os.path.isfile(os.path.join(software, 'final.mtz')):
                os.system('ln -s {0!s} init.mtz'.format(os.path.join(software, 'final.mtz')))
            if os.path.isfile(os.path.join(software, 'final.mmcif')):
                os.system('ln -s {0!s} init.mmcif'.format(os.path.join(software, 'final.mmcif')))
            if os.path.isfile(os.path.join(software, 'reindexed.mtz')):
                os.system('ln -s {0!s} free.mtz'.format(os.path.join(software, 'reindexed.mtz')))

"""

def link_init_refinement_files(logger, projectDir, sample, initpdb, initmtz, initcif, freemtz):
    os.chdir(os.path.join(projectDir, '2-initial_refine', sample))
    if os.path.isfile(initpdb):
        os.system('ln -s {0!s} init.pdb'.format(initpdb))
    if os.path.isfile(initmtz):
        os.system('ln -s {0!s} init.mtz'.format(initmtz))
    if os.path.isfile(initcif):
        os.system('ln -s {0!s} init.mmcif'.format(initcif))
    if os.path.isfile(freemtz):
        os.system('ln -s {0!s} free.mtz'.format(freemtz))


def get_refinement_files(logger, projectDir, sample, software):
    d = {}
    initpdb = None
    initmtz = None
    initcif = None
    freemtz = None
    os.chdir(os.path.join(projectDir, '2-initial_refine', sample))
    if software == 'dimple':
        if os.path.isfile(os.path.join(software, 'final.pdb')):
            initpdb = os.path.relpath(os.path.realpath(os.path.join(software, 'final.pdb')))
        if os.path.isfile(os.path.join(software, 'final.mtz')):
            initmtz = os.path.relpath(os.path.realpath(os.path.join(software, 'final.mtz')))
        if os.path.isfile(os.path.join(software, 'final.mmcif')):
            initcif = os.path.relpath(os.path.realpath(os.path.join(software, 'final.mmcif')))
        if os.path.isfile(os.path.join(software, 'reindexed.mtz')):
            freemtz = os.path.relpath(os.path.realpath(os.path.join(software, 'reindexed.mtz')))
    elif software == 'pipedream':
        if os.path.isfile(os.path.join(software, 'refine', 'refine.pdb')):
            initpdb = os.path.relpath(os.path.realpath(os.path.join(software, 'refine', 'refine.pdb')))
        if os.path.isfile(os.path.join(software, 'refine', 'refine.mtz')):
            initmtz = os.path.relpath(os.path.realpath(os.path.join(software, 'refine', 'refine.mtz')))
        if os.path.isfile(os.path.join(software, 'refine', 'BUSTER_model.cif')):
            initcif = os.path.relpath(os.path.realpath(os.path.join(software, 'refine', 'BUSTER_model.cif')))
        if initpdb:
            for line in open(initpdb):
                if line.startswith("REMARK   MTZ "):
                    freemtz = os.path.relpath(os.path.realpath(line.split()[3]))
    else:
        logger.error('auto-refinement pipeline does not exist')
    d['mounted_crystal_code'] = sample
    d['initial_refinement_pipeline'] = software
    d['initial_refinement_pdb_file'] = initpdb
    d['initial_refinement_mtz_file'] = initmtz
    d['initial_refinement_cif_file'] = initcif
    d['rfree_mtz_file'] = freemtz
    d['initial_refinement_directory'] = os.path.join(projectDir, '2-initial_refine', sample, software)
    return initpdb, initmtz, initcif, freemtz, d


def get_db_dict_from_model_cif(logger, model_cif, d, software):
    logger.info('reading information from {0!s}'.format(model_cif))
    doc = gemmi.cif.read_file(model_cif)
    for block in doc:
        if block.find_pair('_refine.ls_R_factor_R_work'):
            d['refine_ls_R_factor_R_work'] = block.find_pair('_refine.ls_R_factor_R_work')[1]
        if block.find_pair('_refine.ls_R_factor_R_free'):
            d['refine_ls_R_factor_R_free'] = block.find_pair('_refine.ls_R_factor_R_free')[1]
        if block.find_pair('_refine.ls_d_res_high'):
            d['refine_ls_d_res_high'] = block.find_pair('_refine.ls_d_res_high')[1]
        if block.find_pair('_refine.ls_d_res_low'):
            d['refine_ls_d_res_low'] = block.find_pair('_refine.ls_d_res_low')[1]
        if block.find_pair('_cell.length_a'):
            d['cell_length_a'] = block.find_pair('_cell.length_a')[1]
        if block.find_pair('_cell.length_b'):
            d['cell_length_b'] = block.find_pair('_cell.length_b')[1]
        if block.find_pair('_cell.length_c'):
            d['cell_length_c'] = block.find_pair('_cell.length_c')[1]
        if block.find_pair('_cell.angle_alpha'):
            d['cell_angle_alpha'] = block.find_pair('_cell.angle_alpha')[1]
        if block.find_pair('_cell.angle_beta'):
            d['cell_angle_beta'] = block.find_pair('_cell.angle_beta')[1]
        if block.find_pair('_cell.angle_gamma'):
            d['cell_angle_gamma'] = block.find_pair('_cell.angle_gamma')[1]
        if block.find_pair('_symmetry.space_group_name_H-M'):
            d['sym_space_group'] = block.find_pair('_symmetry.space_group_name_H-M')[1]
        if block.find_pair('_symmetry.Int_Tables_number'):
            d['sym_Int_Tables_number'] = block.find_pair('_symmetry.Int_Tables_number')[1]
        if block.find_loop('_refine_ls_restr.type'):
            if software == 'dimple':
                table = block.find('_refine_ls_restr.', ['type', 'dev_ideal'])
                d['refine_r_bond_refined_d'] = str(list(table.find_row('r_bond_refined_d'))[1])
                d['refine_r_angle_refined_deg'] = str(list(table.find_row('r_angle_refined_deg'))[1])
            if software == 'pipedream':
                table = block.find('_refine_ls_restr.', ['type', 'dev_ideal'])
                d['refine_r_bond_refined_d'] = str(list(table.find_row('t_bond_d'))[1])
                d['refine_r_angle_refined_deg'] = str(list(table.find_row('t_angle_deg'))[1])
        if 'refine_ls_R_factor_R_free' in d:
            if float(d['refine_ls_R_factor_R_free']) < 0.3:
                d['initial_refinement_outcome'] = "success"
            elif float(d['refine_ls_R_factor_R_free']) >= 0.3 and float(d['refine_ls_R_factor_R_free']) < 0.35:
                d['initial_refinement_outcome'] = "success - elevated Rfree"
            else:
                d['initial_refinement_outcome'] = "fail - Rfree too high"
        return d
