import os
import time
import gemmi
from PIL import Image
import json
import logging
import glob
import sys
from bz2 import BZ2File as bzopen


def init_logger(logfile):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s | %(levelname)s | %(message)s', '%m-%d-%Y %H:%M:%S')

    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.DEBUG)
    stdout_handler.setFormatter(formatter)
    logger.addHandler(stdout_handler)

    file_handler = logging.FileHandler(logfile)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


def get_sample_list(logger, fragmaxcsv):
    sampleList = []
    for line in open(fragmaxcsv):
        sample = line.split(',')[0]
        sampleList.append(sample)
    logger.info('found {0!s} samples in summary csv file\n'.format(len(sampleList)))
    return sampleList


def get_proposal_and_session_and_protein(processDir):
    proposal = processDir.split('/')[4]
    session = processDir.split('/')[5]
    protein = processDir.split('/')[7]
    return proposal, session, protein


def create_sample_folder(logger, projectDir, sample):
    os.chdir(os.path.join(projectDir, '1-process'))
    if not os.path.isdir(sample):
        logger.info('creating folder {0!s} in 1-process')
        os.mkdir(sample)


def prepare_folders_and_files(logger, projectDir, sample, proposal, session, run, protein, processDir):
    create_proposal_session_run_folder(logger, projectDir, sample, proposal, session, run)
    create_image_folder(logger, projectDir, sample, proposal, session, run)
    find_crystal_snapshots(logger, projectDir, sample, proposal, session, protein, run)
    find_dozor_plot(logger, processDir, projectDir, sample, proposal, session, run)


def create_proposal_session_run_folder(logger, projectDir, sample, proposal, session, run):
    os.chdir(os.path.join(projectDir, '1-process', sample))
    if not os.path.isdir('{0!s}-{1!s}-{2!s}'.format(proposal, session, run)):
        logger.info('creating folder {0!s}-{1!s}-{2!s} in 1-process/{3!s}'.format(proposal, sample, run, sample))
        os.mkdir('{0!s}-{1!s}-{2!s}'.format(proposal, session, run))


def create_pipeline_folder(logger, projectDir, sample, proposal, session, run, pipeline):
    os.chdir(os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run)))
    if not os.path.isdir(pipeline):
        logger.info('creating folder {0!s}'.format(pipeline))
        os.mkdir(pipeline)


def create_image_folder(logger, projectDir, sample, proposal, session, run):
    os.chdir(os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run)))
    if not os.path.isdir('images'):
        logger.info('creating images folder')
        os.mkdir('images')


def find_crystal_snapshots(logger, projectDir, sample, proposal, session, protein, run):
    logger.info('looking for crystal snapshots...')
    os.chdir(os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run), 'images'))
    for img in glob.glob('/data/staff/ispybstorage/pyarch/visitors/{0!s}/{1!s}/raw/{2!s}/{3!s}/*.snapshot.jpeg'.format(
            proposal, session, protein, sample)):
        logger.info('copying {0!s}'.format(img))
        os.system('/bin/cp {0!s} .'.format(img))


def find_dozor_plot(logger, processDir, projectDir, sample, proposal, session, run):
    logger.info('looking for dozor plots...')
    os.chdir(os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run), 'images'))
    for img in glob.glob(os.path.join(processDir, sample, run, 'ImgQIndicator_proc', 'cn*', 'ControlPyDozor*', 'dozor_*.png')):
        logger.info('copying {0!s}'.format(img))
        os.system('/bin/cp {0!s} dozor.png'.format(img))


def get_processing_pipelines():
    pipelines = [
        'autoproc',
        'staraniso',
        'xia2dials'
#        'xia2xds'
    ]
    return pipelines


def get_pipeline_path(pipeline):
    mtzpath = None
    mtz_extension = None
    log_extension = None
    cif_extension = None
    if pipeline == 'autoproc':
        mtzpath = os.path.join('autoPROC', 'cn*', 'AutoPROCv1_*noanom*', 'HDF5_1', 'truncate-unique.mtz')
        mtz_extension = 'HDF5_1/truncate-unique.mtz'
        log_extension = 'HDF5_1/aimless.log'
        cif_extension = 'Data_2_autoPROC_TRUNCATE_all.cif'
    elif pipeline == 'staraniso':
        mtzpath = os.path.join('autoPROC', 'cn*', 'AutoPROCv1_*noanom*', 'HDF5_1', 'staraniso_alldata-unique.mtz')
        mtz_extension = 'HDF5_1/staraniso_alldata-unique.mtz'
        log_extension = 'HDF5_1/staraniso_alldata.log'
        cif_extension = 'Data_1_autoPROC_STARANISO_all.cif'
    elif pipeline == 'xia2dials':
        mtzpath = os.path.join('xia2DIALS', 'cn*', 'Xia2DIALSv1_*noanom', 'DataFiles', 'AUTOMATIC_DEFAULT_free.mtz')
        mtz_extension = 'DataFiles/AUTOMATIC_DEFAULT_free.mtz'
        log_extension = 'LogFiles/AUTOMATIC_DEFAULT_SCALE.log'
        cif_extension = 'DataFiles/xia2.mmcif.bz2'
    return mtzpath, mtz_extension, log_extension, cif_extension


def get_process_files(logger, mtzfile, projectDir, sample, proposal, session,
                      run, pipeline, collection_date, mtz_extension, cif_extension, log_extension, status):
    logfile = None
    ciffile = None
    mtz = mtz_info(mtzfile)
    wavelength = mtz['wavelength']
    create_pipeline_folder(logger, projectDir, sample, proposal, session, run, pipeline)
    if os.path.isfile(mtzfile.replace(mtz_extension, log_extension)):
        logger.info('found LOG file')
        logfile = mtzfile.replace(mtz_extension, log_extension)
    else:
        logger.error('cannot find LOG file')
    if os.path.isfile(mtzfile.replace(mtz_extension, cif_extension)):
        logger.info('found CIF file')
        ciffile = mtzfile.replace(mtz_extension, cif_extension)
    else:
        logger.error('cannot find CIF file')
    if logfile and ciffile:
        copy_files_to_project_folder(projectDir, sample, run, proposal, session, pipeline,
                                                mtzfile, logfile, ciffile, collection_date, wavelength)
    else:
        logger.error('MTZ file exists, but either LOG or CIF file missing')
    status = get_status(mtzfile, status)
    return status


def get_timestamp_from_master_file(sample_folder, run):
    raw_folder = sample_folder.replace("process", "raw")
    master_file = run[4:-1] + "master.h5"
    master = os.path.join(raw_folder, master_file)
    t_c = os.path.getctime(master)
    c_t = time.ctime(t_c)
    t_obj = time.strptime(c_t)
    t_stamp = time.strftime("%Y-%m-%d", t_obj)
    return t_stamp, master


def add_biomax_mmcif_header_items(wavelength, collection_date):
    # date: e.g. 2019-08-02
    header = (
        "data_process_cif\n"
        "\n"
        "_diffrn.id 1\n"
        "_diffrn.crystal_id 1\n"
        "_diffrn.ambient_temp 100\n"
        "\n"
        "loop_\n"
        "_diffrn_source.source\n"
        "_diffrn_source.diffrn_id\n"
        "_diffrn_source.pdbx_synchrotron_site\n"
        "_diffrn_source.type\n"
        "_diffrn_source.pdbx_wavelength_list\n"
        "SYNCHROTRON 1 'MAX IV' 'MAX IV BEAMLINE BioMAX' {0!s}\n".format(wavelength) +
        "\n"
        "loop_\n"
        "_diffrn_detector.type\n"
        "_diffrn_detector.detector\n"
        "_diffrn_detector.pdbx_collection_date\n"
        "_diffrn_detector.diffrn_id\n"
        "'DECTRIS EIGER X 16M' PIXEL {0!s} 1\n".format(collection_date) +
        "\n"
        "loop_\n"
        "_diffrn_radiation.wavelength_id\n"
        "_diffrn_radiation.pdbx_scattering_type\n"
        "_diffrn_radiation.pdbx_diffrn_protocol\n"
        "_diffrn_radiation.diffrn_id\n"
        "_diffrn_radiation.pdbx_monochromatic_or_laue_m_l\n"
        "1 x-ray 'SINGLE WAVELENGTH' 1 M\n"
    )
    return header


def write_mmcif_header(cif, cif_name, collection_date, wavelength):
    cifLines = add_biomax_mmcif_header_items(wavelength, collection_date)
    previous_line = ''
    if cif_name.endswith('.bz2'):
        for line in bzopen(cif):
            if '_pdbx_diffrn_unmerged_cell' in line:
                break
            else:
                cifLines += previous_line
            previous_line = line
    else:
        for line in open(cif):
            if line.startswith('_refln.'):
                break
            else:
                cifLines += previous_line
            previous_line = line
    f = open(cif_name.replace(' ', '').replace('.bz2', ''), 'w')
    f.write(cifLines)
    f.close()


def copy_files_to_project_folder(projectDir, sample, run, proposal, session, pipeline,
                                 mtz, log, cif, collection_date, wavelength):
    os.chdir(os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run), pipeline))
    mtz_name = mtz.split('/')[len(mtz.split('/'))-1]
    log_name = log.split('/')[len(log.split('/'))-1]
    cif_name = cif.split('/')[len(cif.split('/'))-1]
    if not os.path.isfile(mtz_name):
        os.system('/bin/cp {0!s} .'.format(mtz))
    if not os.path.isfile(log_name):
        os.system('/bin/cp {0!s} .'.format(log))
    if not os.path.isfile(cif_name):
        write_mmcif_header(cif, cif_name, collection_date, wavelength)
    create_process_symlink(mtz_name, log_name, cif_name)


def create_process_symlink(mtz_name, log_name, cif_name):
    if not os.path.isfile('process.mtz'):
        os.system('ln -s {0!s} process.mtz'.format(mtz_name))
    if not os.path.isfile('process.log'):
        os.system('ln -s {0!s} process.log'.format(log_name))
    if not os.path.isfile('process.cif'):
        os.system('ln -s {0!s} process.cif'.format(cif_name.replace('.bz2', '')))



def mtz_info(mtzfile):
    mtzDict = {}
    mtz = gemmi.read_mtz_file(mtzfile)
    mtzDict['unitcell'] = mtz.cell
    mtzDict['unitcell_volume'] = mtz.cell.volume
    mtzDict['point_group'] = mtz.spacegroup.point_group_hm()
    mtzDict['resolution_high'] = mtz.resolution_high()
    mtzDict['space_group'] = mtz.spacegroup.hm
    mtzDict['wavelength'] = mtz.dataset(1).wavelength
#    mtzDict['lattice'] = mtz.spacegroup.crystal_system_str()
    mtzDict['lattice'] = mtz.spacegroup.hm[0]
    return mtzDict


def cif_info(ciffile):
    cifDict = {}
    doc = gemmi.cif.read_file(ciffile)
    for block in doc:
        if block.find_pair('_symmetry.space_group_name_H-M'):
            cifDict['space_group'] = str(block.find_pair('_symmetry.space_group_name_H-M')[1])
        if block.find_pair('_cell.length_a'):
            cifDict['a'] = str(block.find_pair('_cell.length_a')[1])
            cifDict['b'] = str(block.find_pair('_cell.length_b')[1])
            cifDict['c'] = str(block.find_pair('_cell.length_c')[1])
            cifDict['alpha'] = str(block.find_pair('_cell.angle_alpha')[1])
            cifDict['beta'] = str(block.find_pair('_cell.angle_beta')[1])
            cifDict['gamma'] = str(block.find_pair('_cell.angle_gamma')[1])
        if block.find_pair('_reflns.d_resolution_low'):
            cifDict['reso_low'] = str(round(float(block.find_pair('_reflns.d_resolution_low')[1]), 2))
#            cifDict['reso_low'] = str(block.find_pair('_reflns.d_resolution_low')[1])
            cifDict['reso_high'] = str(round(float(block.find_pair('_reflns.d_resolution_high')[1]), 2))
#        if block.find_pair('_reflns.pdbx_netI_over_sigmaI'):
            cifDict['pdbx_netI_over_sigmaI'] = str(block.find_pair('_reflns.pdbx_netI_over_sigmaI')[1])
            cifDict['pdbx_redundancy'] = str(block.find_pair('_reflns.pdbx_redundancy')[1])
            if 'staraniso' in ciffile:
                cifDict['percent_possible_obs'] = str(block.find_pair('_reflns.pdbx_percent_possible_spherical')[1])
            else:
                cifDict['percent_possible_obs'] = str(block.find_pair('_reflns.percent_possible_obs')[1])
            cifDict['pdbx_number_measured_all'] = str(block.find_pair('_reflns.pdbx_number_measured_all')[1])
            cifDict['pdbx_CC_half'] = str(block.find_pair('_reflns.pdbx_CC_half')[1])

        if block.find_loop('_reflns_shell.Rmerge_I_obs'):
            Rmerge_I_obs = list(block.find_loop('_reflns_shell.Rmerge_I_obs'))
            cifDict['Rmerge_I_obs_low'] = Rmerge_I_obs[0]

        if block.find_loop('_reflns_shell.meanI_over_sigI_obs'):
            meanI_over_sigI_obs = list(block.find_loop('_reflns_shell.meanI_over_sigI_obs'))
            cifDict['meanI_over_sigI_obs_high'] = meanI_over_sigI_obs[len(meanI_over_sigI_obs)-1]
    return cifDict


def get_status(mtzfile, status):
    if not status is 'OK':
        if mtzfile:
            mtzDict = mtz_info(mtzfile)
            if float(mtzDict['resolution_high']) < 2.5:
                status = 'OK'
            elif float(mtzDict['resolution_high']) >= 2.8 and float(mtzDict['resolution_high']) < 3.2:
                if status == 'FAIL - low resolution' or status == 'FAIL - no processing result':
                    status = 'OK - medium resolution'
            else:
                if status == 'FAIL - no processing result':
                    status = 'FAIL - low resolution'
        else:
            status = 'FAIL - no processing result'
    return status


def make_thumbnail(folder, image):
    # save a smaller version of the image
    thumbnail = image.replace('.png', '_thumb.png')
#    width_30 = int(round(width_100 * 0.3, 0))
    img = Image.open(os.path.join(folder, image))
#    wpercent = (width_30/float(width_100))
#    hsize = int((float(height_100)*float(wpercent)))
#    img = img.resize((width_30,hsize), Image.ANTIALIAS)
    img.save(os.path.join(folder, thumbnail))


def write_json_info_file(logger, projectDir, sample, collection_date, run, proposal, session,
                         protein, status, master, pipeline):
    os.chdir(os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run)))
    if pipeline:
        os.chdir(pipeline)
    d = {}
    d['sample'] = sample
    d['collection_date'] = collection_date
    d['run'] = run
    d['proposal'] = proposal
    d['session'] = session
    d['protein'] = protein
    d['status'] = status
    d['master'] = master
    d['pipeline'] = pipeline
    logger.info('saving info.json file...')
    with open("info.json", "w") as outfile:
        json.dump(d, outfile)


def read_reference_pdb_files(logger, projectDir):
    logger.info('reading point group, lattice and unit cell volume from reference pdb files...')
    os.chdir(os.path.join(projectDir, '0-model'))
    ref_dict = {}
    for pdbfile in glob.glob('*.pdb'):
        structure = gemmi.read_pdb(pdbfile)
        unitcell = structure.cell
        unitcell_volume = unitcell.volume
        sym = gemmi.find_spacegroup_by_name(structure.spacegroup_hm)
        point_group = sym.point_group_hm()
        lattice = sym.centring_type()
        ref_dict[pdbfile] = [point_group, lattice, unitcell_volume]
        logger.info('{0!s} - pg: {1!s}, lattice: {2!s}, uc_vol: {3!s}'.format(
            pdbfile, point_group, lattice, unitcell_volume))
    if not ref_dict:
        logger.warning('could not find any PDB files in 0-model')
    return ref_dict


def retain_results_with_similar_ucvol_and_pg_as_ref_pdb(logger, proc_dict, ref_dict):
    logger.info('checking if PDB files in 0-model folder are similar to the processing MTZ files...')
    match_dict = {}
    if ref_dict:
        for f in proc_dict:
            pgr_mtz = proc_dict[f]['point_group']
            ucv_mtz = float(proc_dict[f]['unitcell_volume'])
            lat_mtz = proc_dict[f]['lattice']
            for p in ref_dict:
                pgr_pdb = ref_dict[p][0]
                ucv_pdb = float(ref_dict[p][2])
                lat_pdb = ref_dict[p][1]
                ucv_diff = abs((ucv_mtz-ucv_pdb))/ucv_pdb
                if pgr_mtz == pgr_pdb and lat_mtz == lat_pdb and ucv_diff < 0.1:
                    logger.info('lattice, point group and unit cell volume of MTZ file matches {0!s}'.format(p))
                    match_dict[f] = proc_dict[f]
                    break
    else:
        logger.warning('seems like no reference PDB files were provided; skipping this selection step...')
    if match_dict:
        proc_dict = match_dict
    else:
        logger.warning('none of the PDB files appears to be similar to any of the auto-processing MTZ files')
    return proc_dict


def retain_results_with_good_low_reso_rmerge(logger, proc_dict):
    logger.info('checking of auto-processing MTZ files have low resolution Rmerge values below 10%...')
    match_dict = {}
    for f in proc_dict:
        if float(proc_dict[f]['Rmerge_I_obs_low']) < 0.1:
            logger.info('{0!s} - Rmerge(low): {1!s}'.format(f, proc_dict[f]['Rmerge_I_obs_low']))
            match_dict[f] = proc_dict[f]
        else:
            logger.error('{0!s} - Rmerge(low): {1!s}'.format(f, proc_dict[f]['Rmerge_I_obs_low']))
    if match_dict:
        proc_dict = match_dict
    else:
        logger.warning('did not find any MTZ file with Rmerge (low) below 10%')
        logger.info('will continue anyway...')
    return proc_dict


def retain_results_which_fit_selection_criterion(logger, proc_dict, select_criterion):
    logger.info('selecting auto-processing results based on {0!s}...'.format(select_criterion))
    match_list = []
    backup_list = []
    for f in proc_dict:
        reso_high = proc_dict[f]['reso_high']
        if select_criterion.startswith('reso'):
            logger.info('added {0!s} with high resolution limit of {1!s} A'.format(f, reso_high))
            match_list.append([f, float(reso_high)])
        elif proc_dict[f]['pipeline'] == select_criterion:
            logger.info('added {0!s} with high resolution limit of {1!s} A'.format(f, reso_high))
            match_list.append([f, float(reso_high)])
        else:
            logger.warning('MTZ does not match criteria, but added {0!s} with high resolution limit of {1!s} A'.format(
                f, reso_high))
            backup_list.append([f, float(reso_high)])
    if not match_list:
        logger.warning('none of the MTZ files fulfilled the selection criteria, but will select the one with highest resolution')
        match_list = backup_list
    logger.info('selecting MTZ file based on highest resolution...')
    bestcif = min(match_list, key=lambda x: x[1])[0]
    logger.info('--> {0!s}'.format(bestcif))
    return bestcif


def link_process_results(logger, projectDir, sample, bestcif):
    logger.info('creating symlinks in {0!s}'.format(os.path.join(projectDir, "1-process", sample)))
    os.chdir(os.path.join(projectDir, "1-process", sample))
    pipeline = bestcif.split('/')[1]
#    json_info = bestcif.replace(pipeline +'/process.cif', 'info.json')
    json_info = bestcif.replace('process.cif', 'info.json')
    if not os.path.isdir('process.mtz'):
        os.system('ln -s {0!s} .'.format(bestcif.replace('.cif', '.mtz')))
    if not os.path.isdir('process.log'):
        os.system('ln -s {0!s} .'.format(bestcif.replace('.cif', '.log')))
    if not os.path.isdir('process.cif'):
        os.system('ln -s {0!s} .'.format(bestcif))
    if not os.path.isdir('info.json'):
        os.system('ln -s {0!s} .'.format(json_info))


def link_info_json_file(logger, projectDir, sample):
    os.chdir(os.path.join(projectDir, "1-process", sample))
    for f in glob.glob(os.path.join('*', 'info.json')):
        logger.info('linking {0!s} to sample directory')
        os.system('ln -s {0!s} .'.format(f))
        

def remove_process_symlinks(logger, projectDir, sample):
    logger.warning('will remove existing process symlinks in {0!s}'.format(os.path.join(
        projectDir, "1-process", sample)))
    os.chdir(os.path.join(projectDir, '1-process', sample))
    os.system('/bin/rm process.*')
    os.system('/bin/rm info.json')


def skip_sample_if_already_selected(logger, projectDir, sample, sample_folder, overwrite):
    process_not_exists = False
    if os.path.isfile(os.path.join(sample_folder, 'process.cif')):
        logger.warning('linked processing results exist, i.e. results for this sample have been selected before')
        if not overwrite:
            logger.info('will skip selection for sample; choose overwrite option if you want to reselect')
            process_not_exists = True
        else:
            logger.warning('you chose to overwrite existing links; will re-select accordingly')
            remove_process_symlinks(logger, projectDir, sample)
    else:
        if overwrite:
            # in case of broken links
            remove_process_symlinks(logger, projectDir, sample)
    return process_not_exists


def read_data_collection_stats(logger, ciffile, proc_dict):
    subfolder = ciffile.split('/')[0]
    pipeline = ciffile.split('/')[1]
    mtzfile = os.path.join(ciffile.replace('.cif', '.mtz'))
    if os.path.isfile(mtzfile):
        logger.info('found {0!s} file'.format(mtzfile))
        mtz = mtz_info(mtzfile)
        mtz['pipeline'] = pipeline
        mtz['subfolder'] = subfolder
        cif = cif_info(ciffile)
        mtz.update(cif)
        proc_dict[os.path.join(ciffile)] = mtz
    else:
        logger.error('process.mtz does not exist in folder')
    return proc_dict


def start_logging(logger, script):
    logger.info('===================================================================================')
    logger.info('>>>>> starting {0!s}...'.format(script))
    logger.info('===================================================================================')


def start_get_autoprocessing_results(logger):
    logger.info('===================================================================================')
    logger.info('>>> START: collating auto-processing results')
    logger.info('===================================================================================')


def end_get_autoprocessing_results(logger):
    logger.info('===================================================================================')
    logger.info('>>> END: finished collating auto-processing results')
    logger.info('===================================================================================')


def start_select_results(logger):
    logger.info('===================================================================================')
    logger.info('>>> START: selecting auto-processing results')
    logger.info('===================================================================================')


def end_select_results(logger):
    logger.info('===================================================================================')
    logger.info('>>> END: finished selecting auto-processing results')
    logger.info('===================================================================================')


def report_parameters(logger, processDir, projectDir, fragmaxcsv, select, select_criterion, overwrite):
    logger.info('> process directory: {0!s}'.format(processDir))
    logger.info('> project directory: {0!s}'.format(projectDir))
    logger.info('> fragmax csv file:  {0!s}'.format(fragmaxcsv))
    logger.info('> select datasets:   {0!s}'.format(select))
    logger.info('> select criterion:  {0!s}'.format(select_criterion))
    logger.info('> overwrite:         {0!s}'.format(overwrite))
    logger.info('===================================================================================')


def check_if_process_directory_exists(logger, processDir, select):
    passed = True
    if not select:
        logger.info('checking process directory: {0!s}'.format(processDir))
        if os.path.isdir(processDir):
            logger.info('process directory exists')
        else:
            logger.error('process directory does not exist')
            passed = False
    return passed


def check_if_project_directory_exists(logger, projectDir, passed):
    logger.info('checking project directory: {0!s}'.format(projectDir))
    if os.path.isdir(projectDir):
        logger.info('project directory exists')
    else:
        logger.error('project directory does not exist')
        passed = False
    return passed


def check_if_fragmaxcsv_exists(logger, fragmaxcsv, passed):
    logger.info('checking fragmaxcsv: {0!s}'.format(fragmaxcsv))
    if os.path.isfile(fragmaxcsv):
        logger.info('fragmaxcsv exists')
    else:
        logger.error('fragmaxcsv does not exist')
        passed = False
    return passed


def check_pdbfiles_in_pdbdir(logger, projectDir, passed):
    logger.info('looking for PDB files in PDB directory...')
    found = False
    for pdb in glob.glob(os.path.join(projectDir, '0-model', '*.pdb')):
        pdbFile = pdb.split('/')[len(pdb.split('/'))-1]
        logger.info('found {0!s}'.format(pdbFile))
        foundCryst = False
        for line in open(pdb):
            if line.startswith('CRYST'):
                logger.info('found CRYST card: {0!s}'.format(line.replace('\n', '')))
                foundCryst = True
                found = True
        if not foundCryst:
            logger.warning('{0!s} does not seem to contain a CRYST card'.format(pdbFile))
    if found:
        logger.info('found at least one PDB file with a CRYST card')
    else:
        logger.warning('did not find a PDB file with a valid CRYST card')
        passed = False
    return passed


def check_select_criterion(logger, select_criterion, passed):
    logger.info('checking selection option: {0!s}'.format(select_criterion))
    supported_options = ['xia2dials', 'xia2xds', 'autoproc', 'staraniso', 'resolution']
    if select_criterion in supported_options:
        logger.info('option exists')
    else:
        logger.error('ERROR: option does not exist')
        passed = False
    return passed


def run_checks(logger, processDir, projectDir, fragmaxcsv, select, select_criterion):
    logger.info('checking input file and command line options...')
    passed = check_if_process_directory_exists(logger, processDir, select)
    passed = check_if_project_directory_exists(logger, projectDir, passed)
    passed = check_if_fragmaxcsv_exists(logger, fragmaxcsv, passed)
    passed = check_pdbfiles_in_pdbdir(logger, projectDir, passed)
    passed = check_select_criterion(logger, select_criterion, passed)
    logger.info('===================================================================================')
    return passed


def usage():
    usage = (
        '\n'
        'usage:\n'
        'ccp4-python 1-process.py -i <process_dir> -o <project_dir> -f <fragmax_csv_file>\n'
        '\n'
        'additional command line options:\n'
        '--input, -i\n'
        '    process directory\n'
        '--output, -o\n'
        '    project directory\n'
        '--fragmaxcsv, -f\n'
        '    fragmax summary csv file\n'
        '--select, -r\n'
        '    select auto-processing results\n'
        '--crtierion, -c\n'
        '    selection criterion (resolution [default], xia2dials, xia2xds, autoproc, staraniso)\n'
        '--overwrite, -o\n'
        '    flag to overwrite selected files\n'
    )
    print(usage)


def check_if_to_continue(logger):
    if sys.version[0] == '2':
        q = raw_input("\n>>> Do you want to continue? (y/n) ")
    else:
        q = input("\n>>> Do you want to continue? (y/n)")
    if not q.lower() == 'y':
        logger.info('you chose not to conutinue at this point; exciting program...')
        sys.exit(2)


