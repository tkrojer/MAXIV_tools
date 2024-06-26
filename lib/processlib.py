import os
import time
import gemmi
from PIL import Image
import json
import logging
import glob
import sys
from bz2 import BZ2File as bzopen
from datetime import datetime

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'lib')))
import processdb


def max_allowed_Rmerge_I_obs_low():
    Rmerge_I_obs_low = 0.25
    return Rmerge_I_obs_low


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
        sample = line.split(',')[0].replace('\n', '')
        sampleList.append(sample)
    logger.info('found {0!s} samples in summary csv file\n'.format(len(sampleList)))
    return sampleList


def get_proposal_and_session_and_protein(processDir):
    try:
        proposal = processDir.split('/')[4]
        session = processDir.split('/')[5]
        protein = processDir.split('/')[7]
        beamline = processDir.split('/')[3]
        category = processDir.split('/')[2]
    except AttributeError:
        proposal = None
        session = None
        protein = None
        beamline = None
        category = None
    return proposal, session, protein, beamline, category

def get_none_proposal_and_session_and_protein():
    proposal = None
    session = None
    protein = None
    beamline = None
    category = None
    return proposal, session, protein, beamline, category

def create_sample_folder(logger, projectDir, sample):
    os.chdir(os.path.join(projectDir, '1-process'))
    if not os.path.isdir(sample):
        logger.info('creating folder {0!s} in 1-process'.format(sample))
        os.mkdir(sample)


def prepare_folders_and_files(logger, projectDir, sample, proposal, session, run, protein, processDir, category, beamline):
    create_proposal_session_run_folder(logger, projectDir, sample, proposal, session, run)
    create_image_folder(logger, projectDir, sample, proposal, session, run)
    crystal_snapshot_list = find_crystal_snapshots(logger, projectDir, sample, proposal, session, protein, run, category, beamline)
    dozor_plot = find_dozor_plot(logger, processDir, projectDir, sample, proposal, session, run)
    return dozor_plot, crystal_snapshot_list

def create_proposal_session_run_folder(logger, projectDir, sample, proposal, session, run):
    os.chdir(os.path.join(projectDir, '1-process', sample))
    if not os.path.isdir('{0!s}-{1!s}-{2!s}'.format(proposal, session, run)):
        logger.info('creating folder {0!s}-{1!s}-{2!s} in 1-process/{3!s}'.format(proposal, sample, run, sample))
        os.mkdir('{0!s}-{1!s}-{2!s}'.format(proposal, session, run))


def create_pipeline_folder(logger, projectDir, sample, proposal, session, run, pipeline, overwrite):
    cont = True
    os.chdir(os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run)))
    if not os.path.isdir(pipeline):
        logger.info('creating folder {0!s}'.format(pipeline))
        os.mkdir(pipeline)
    else:
        logger.warning('folder {0!s} exists'.format(pipeline))
        if overwrite:
            logger.warning('overwrite is True; removing existing folder and creating a new one')
            os.system('rm -fr {0!s}'.format(pipeline))
            os.mkdir(pipeline)
        else:
            logger.warning('overwrite is False; skipping')
            cont = False
    return cont

def create_image_folder(logger, projectDir, sample, proposal, session, run):
    os.chdir(os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run)))
    if not os.path.isdir('images'):
        logger.info('creating images folder')
        os.mkdir('images')


def find_crystal_snapshots(logger, projectDir, sample, proposal, session, protein, run, category, beamline):
    crystal_snapshot_list = []
    snapshot_dir = os.path.join('/data', 'staff', 'ispybstorage', category, beamline, category, proposal, session, 'raw',
                                protein, sample)
    logger.info('looking for crystal snapshots in {0!s}'.format(snapshot_dir))
    image_dir = os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run), 'images')
    os.chdir(image_dir)

    #/data/staff/ispybstorage/proprietary/biomax/proprietary/20230893/20230624/raw/GEN2110_A/GEN2110_A-x0112/GEN2110_A-x0112_1_1.snapshot.jpeg

#    for img in glob.glob('/data/staff/ispybstorage/pyarch/visitors/{0!s}/{1!s}/raw/{2!s}/{3!s}/*.snapshot.jpeg'.format(
#            proposal, session, protein, sample)):
    logger.info(f"glob string: {run.replace('xds_', '')[:-2]}*.snapshot.jpeg")
    for img in sorted(glob.glob(os.path.join(snapshot_dir, f"{run.replace('xds_', '')[:-2]}*.snapshot.jpeg"))):
        logger.info('copying {0!s}'.format(img))
        os.system('/bin/cp {0!s} .'.format(img))
        img_name = img[img.rfind('/') + 1:]
        if len(crystal_snapshot_list) <= 4:
            crystal_snapshot_list.append(os.path.join(image_dir, img_name))

    if not crystal_snapshot_list:
        logger.warning(f"could not find crystal snapshots for run {run}. This is most likely because this is either a line scan or an automitcally processed dataset")
        logger.info(f"trying to go backwards in run number to find runs from optical prealignment")
        current_run = run.replace('xds_', '')[:-2]
        run_number = int(current_run.split('_')[len(current_run.split('_'))-1])
        run_base = current_run[:current_run.rfind('_')]
        for n in reversed(range(run_number)):
            for img in sorted(glob.glob(os.path.join(snapshot_dir, f"{run_base}_{n}*.snapshot.jpeg"))):
                logger.info('copying {0!s}'.format(img))
                os.system('/bin/cp {0!s} .'.format(img))
                img_name = img[img.rfind('/') + 1:]
                if len(crystal_snapshot_list) <= 4:
                    crystal_snapshot_list.append(os.path.join(image_dir, img_name))
            if crystal_snapshot_list:
                break
    while len(crystal_snapshot_list) < 4:
        crystal_snapshot_list.append('None')    # None is deliberately a string
    return crystal_snapshot_list


def find_dozor_plot(logger, processDir, projectDir, sample, proposal, session, run):
    logger.info('looking for dozor plots...')
    dozor_plot = None
    os.chdir(os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run), 'images'))
    for img in glob.glob(os.path.join(processDir, sample, run, 'ImgQIndicator_proc', 'cn*', 'ControlPyDozor*', 'dozor_*.png')):
        logger.info('copying {0!s}'.format(img))
        os.system('/bin/cp {0!s} dozor.png'.format(img))
        dozor_plot = os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run), 'images', 'dozor.png')
    return dozor_plot


def get_processing_pipelines():
    pipelines = [
        'autoproc',
        'autoproc_old',
        'staraniso',
        'staraniso_old',
        'xia2dials',
        'xia2dials_old',
        'xia2xds'
    ]
    return pipelines

def get_manual_processing_pipelines():
    pipelines = [
        'autoproc_manual'
    ]
    return pipelines

def get_manual_pipeline_name(logger, pipeline, mtzfile):
    manual_pipeline = None
    tmp = mtzfile.split('/')
    for item in tmp:
        if item.startswith(pipeline + '_'):
            manual_pipeline = item
            logger.info('manual processing pipeline: ' + manual_pipeline)
            break
    if not manual_pipeline:
        logger.error('cannot identify name of manual processing pipeline in file ' + mtzfile)
    return manual_pipeline


def get_pipeline_path(pipeline):
    mtzpath = None
    mtz_extension = None
    log_extension = None
    cif_extension = None
    mtz_unmerged = None
    if pipeline == 'autoproc':
        mtzpath = os.path.join('MAXIVFastProcessingTask_0', 'AutoPROCTask_*', 'AutoPROCExecTask_*', 'AutoPROCExec_*', 'truncate-unique.mtz')
        mtz_extension = 'truncate-unique.mtz'
        log_extension = 'aimless.log'
        cif_extension = 'Data_2_autoPROC_TRUNCATE_all.cif'
        mtz_unmerged = 'aimless_unmerged.mtz'
    elif pipeline == 'autoproc_old':
        mtzpath = os.path.join('autoPROC', 'cn*', 'AutoPROCv1_*anom*', 'HDF5_1', 'truncate-unique.mtz')
#        mtzpath = os.path.join('autoPROC', 'cn*', 'AutoPROCv1_*noanom*', 'HDF5_1', 'truncate-unique.mtz')
        mtz_extension = 'HDF5_1/truncate-unique.mtz'
        log_extension = 'HDF5_1/aimless.log'
        cif_extension = 'Data_2_autoPROC_TRUNCATE_all.cif'
        mtz_unmerged = 'HDF5_1/aimless_unmerged.mtz'
    elif pipeline == 'autoproc_manual':
    #        mtzpath = os.path.join('autoPROC', 'cn*', 'AutoPROCv1_*noanom*', 'HDF5_1', 'truncate-unique.mtz')
        mtzpath = os.path.join('autoproc_manual', 'HDF5_1', 'truncate-unique.mtz')
        mtz_extension = 'HDF5_1/truncate-unique.mtz'
#        mtzpath = os.path.join('autoproc_manual', 'HDF5_1', 'truncate.mtz')
#        mtz_extension = 'HDF5_1/truncate.mtz'
        log_extension = 'HDF5_1/aimless.log'
        cif_extension = 'Data_2_autoPROC_TRUNCATE_all.cif'
#        cif_extension = 'HDF5_1/aimless.mrfana20.cif'
        mtz_unmerged = 'HDF5_1/aimless_unmerged.mtz'
    elif pipeline == 'staraniso':
        mtzpath = os.path.join('MAXIVFastProcessingTask_0', 'AutoPROCTask_*', 'AutoPROCExecTask_*', 'AutoPROCExec_*', 'staraniso_alldata-unique.mtz')
        mtz_extension = 'staraniso_alldata-unique.mtz'
        log_extension = 'staraniso_alldata.log'
        cif_extension = 'Data_1_autoPROC_STARANISO_all.cif'
        mtz_unmerged = 'aimless_unmerged.mtz'
    elif pipeline == 'staraniso_old':
        mtzpath = os.path.join('autoPROC', 'cn*', 'AutoPROCv1_*anom*', 'HDF5_1', 'staraniso_alldata-unique.mtz')
#        mtzpath = os.path.join('autoPROC', 'cn*', 'AutoPROCv1_*noanom*', 'HDF5_1', 'staraniso_alldata-unique.mtz')
        mtz_extension = 'HDF5_1/staraniso_alldata-unique.mtz'
        log_extension = 'HDF5_1/staraniso_alldata.log'
        cif_extension = 'Data_1_autoPROC_STARANISO_all.cif'
        mtz_unmerged = 'HDF5_1/aimless_unmerged.mtz'
    elif pipeline == 'xia2dials':
##        mtzpath = os.path.join('xia2DIALS', 'cn*', 'Xia2DIALSv1_*noanom', 'DataFiles', 'AUTOMATIC_DEFAULT_free.mtz')
#        mtzpath = os.path.join('xia2DIALS', 'cn*', 'Xia2DIALSv1_*anom', 'DataFiles', 'AUTOMATIC_DEFAULT_free.mtz')
# MAXIVFastProcessingTask_0/Xia2DIALSTask_0/Xia2DialsExecTask_0
        mtzpath = os.path.join('MAXIVFastProcessingTask_0', 'Xia2DIALSTask_*', 'Xia2DialsExecTask_*', 'DataFiles', 'AUTOMATIC_DEFAULT_free.mtz')
        mtz_extension = 'DataFiles/AUTOMATIC_DEFAULT_free.mtz'
        log_extension = 'LogFiles/AUTOMATIC_DEFAULT_SCALE.log'
        cif_extension = 'DataFiles/xia2.mmcif.bz2'
        mtz_unmerged = 'DataFiles/AUTOMATIC_DEFAULT_scaled_unmerged.mtz'
    elif pipeline == 'xia2dials_old':
        mtzpath = os.path.join('xia2DIALS', 'cn*', 'Xia2DIALSv1_*anom', 'DataFiles', 'AUTOMATIC_DEFAULT_free.mtz')
        mtz_extension = 'DataFiles/AUTOMATIC_DEFAULT_free.mtz'
        log_extension = 'LogFiles/AUTOMATIC_DEFAULT_SCALE.log'
        cif_extension = 'DataFiles/xia2.mmcif.bz2'
        mtz_unmerged = 'DataFiles/AUTOMATIC_DEFAULT_scaled_unmerged.mtz'
    elif pipeline == 'xia2xds':
#        mtzpath = os.path.join('xia2DIALS', 'cn*', 'Xia2DIALSv1_*noanom', 'DataFiles', 'AUTOMATIC_DEFAULT_free.mtz')
        mtzpath = os.path.join('xia2XDS', 'cn*', 'Xia2DIALSv1_*anom', 'DataFiles', 'AUTOMATIC_DEFAULT_free.mtz')
        mtz_extension = 'DataFiles/AUTOMATIC_DEFAULT_free.mtz'
        log_extension = 'LogFiles/AUTOMATIC_DEFAULT_SCALE.log'
        cif_extension = 'DataFiles/xia2.mmcif.bz2'
        mtz_unmerged = 'DataFiles/AUTOMATIC_DEFAULT_scaled_unmerged.mtz'
    return mtzpath, mtz_extension, log_extension, cif_extension, mtz_unmerged


def process_files_for_run_pipeline_exist(logger, projectDir, sample, proposal, session, run, pipeline):
    skip = False
    process_cif = os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run),
                               pipeline, 'process.cif')
    if os.path.isfile(process_cif):
        logger.warning('process.cif for {0!s} exists in {1!s}-{2!s}-{3!s}; skipping...'.format(
            pipeline, proposal, session, run))
        skip = True
    return skip


def get_process_files(logger, mtzfile, projectDir, sample, proposal, session,
                      run, pipeline, collection_date, mtz_extension, cif_extension, log_extension, status, mtz_unmerged):
    logfile = None
    ciffile = None
    unm_mtz = None
    mrfana_ciffile = None
    proc_header = mtzfile.replace(mtz_extension, 'process_header.cif')
    mtz = mtz_info(mtzfile)
    wavelength = mtz['wavelength']
    overwrite = False
    create_pipeline_folder(logger, projectDir, sample, proposal, session, run, pipeline, overwrite)

    logger.info('log_extension: {0!s}'.format(log_extension))
    logger.info('mtz_extension: {0!s}'.format(mtz_extension))
    logger.info('mtzfile: {0!s}'.format(mtzfile))

    logger.info('looking for logfile: {0!s}'.format(mtzfile.replace(mtz_extension, log_extension)))
    if os.path.isfile(mtzfile.replace(mtz_extension, log_extension)):
        logfile = mtzfile.replace(mtz_extension, log_extension)
        logger.info('found LOG file: ' + logfile)
    else:
        logger.error('cannot find LOG file')

    logger.info('looking for ciffile' + mtzfile.replace(mtz_extension, cif_extension))
    if os.path.isfile(mtzfile.replace(mtz_extension, cif_extension)):
        ciffile = mtzfile.replace(mtz_extension, cif_extension)
        logger.info('found CIF file: ' + ciffile)
    else:
        logger.error('cannot find CIF file')
        logger.info('seems like something went wrong with autoproc; this may not be a problem but check messages below')

    logger.info('checking HDF5_1 folder for mrfana cif file...')
    logger.info('looking for mrfana ciffile' + str(mtzfile.replace(mtz_extension, "HDF5_1/aimless.mrfana20.cif")))
    if os.path.isfile(mtzfile.replace(mtz_extension, "HDF5_1/aimless.mrfana20.cif")):
        mrfana_ciffile = mtzfile.replace(mtz_extension, "HDF5_1/aimless.mrfana20.cif")
        logger.info('found ' + mrfana_ciffile)
    else:
        logger.error('cannot find ' + mrfana_ciffile)

    logger.info('looking for unmerged mtz file: {0!s}'.format(mtzfile.replace(mtz_extension, mtz_unmerged)))
    if os.path.isfile(mtzfile.replace(mtz_extension, mtz_unmerged)):
        unm_mtz = mtzfile.replace(mtz_extension, mtz_unmerged)
        logger.info('found UNMERGED MTZ file: ' + mtzfile.replace(mtz_extension, mtz_unmerged))
    else:
        logger.error('cannot find UNMERGED MTZ file')

    if logfile and ciffile:
        mtzfile, logfile, ciffile, mrfana_ciffile = copy_files_to_project_folder(logger, projectDir, sample, run, proposal, session, pipeline,
                                                mtzfile, logfile, ciffile, collection_date, wavelength, unm_mtz, mrfana_ciffile)
    else:
        logger.error('MTZ file exists, but either LOG or CIF file missing')

    status = get_status(logger, mtzfile, mtz, ciffile, status, mrfana_ciffile, proc_header)
    logger.info('current status: ' + status)
    return status, logfile, ciffile, mtzfile, mrfana_ciffile


def get_timestamp_from_master_file(logger, sample_folder, run):
    master = None
    t_stamp = None
    create_date = None
    raw_folder = sample_folder.replace("process", "raw")
    master_file = run[4:-1] + "master.h5"
    logger.info('looking for {0!s}'.format(os.path.join(raw_folder, master_file)))
    if os.path.isfile(os.path.join(raw_folder, master_file)):
        logger.info('found .h5 master file')
        master = os.path.join(raw_folder, master_file)
        t_c = os.path.getctime(master)
        # need create_date for database because t_stamp has only date for mmcif file
        create_date = datetime.fromtimestamp(t_c)
        c_t = time.ctime(t_c)
        t_obj = time.strptime(c_t)
        t_stamp = time.strftime("%Y-%m-%d", t_obj)
    else:
        logger.error('cannot find .h5 master file...')
    return t_stamp, master, create_date


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


def write_mmcif_header(logger, cif, cif_name, collection_date, wavelength):
    logger.info('writing mmcif header file...')
    cifLines = add_biomax_mmcif_header_items(wavelength, collection_date)
    previous_line = ''
    if cif_name.endswith('.bz2'):
        for line in bzopen(cif):
#            print(line.decode('ASCII'))
            if '_pdbx_diffrn_unmerged_cell' in str(line):
                break
            else:
                try:
                    cifLines += previous_line.decode('ASCII')
                except AttributeError:
                    cifLines += previous_line
            previous_line = line
#        sys.exit()
    else:
        for line in open(cif):
            if line.startswith('_refln.'):
                break
            else:
                cifLines += previous_line
            previous_line = line
    f = open('process_header.cif', 'w')
    f.write(cifLines)
    f.close()
    logger.warning('done writing mmcif header file')



def copy_files_to_project_folder(logger, projectDir, sample, run, proposal, session, pipeline,
                                 mtz, log, cif, collection_date, wavelength, unm_mtz, mrfana_ciffile):
    os.chdir(os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run), pipeline))
    mtz_name = mtz.split('/')[len(mtz.split('/'))-1]
    log_name = log.split('/')[len(log.split('/'))-1]
    cif_name = cif.split('/')[len(cif.split('/'))-1]
    unm_name = unm_mtz.split('/')[len(unm_mtz.split('/')) - 1]
    mrf_name = mrfana_ciffile.split('/')[len(unm_mtz.split('/')) - 1]
    if not os.path.isfile(mtz_name):
        logger.info(f'copying {mtz}')
        os.system('/bin/cp {0!s} .'.format(mtz))
    if not os.path.isfile(log_name):
        logger.info(f'copying {log}')
        os.system('/bin/cp {0!s} .'.format(log))
    if not os.path.isfile(cif_name):
        logger.info(f'copying {cif}')
        os.system('/bin/cp {0!s} .'.format(cif))
    if not os.path.isfile(mrf_name):
        logger.info(f'copying {mrfana_ciffile}')
        os.system('/bin/cp {0!s} .'.format(mrfana_ciffile))
        check_first_line_of_aimless_mrfana_cif(logger, mrf_name, projectDir, sample, proposal, session, run, pipeline)

    logger.info('unmerged_mtz: {0!s}'.format(unm_mtz))
    if unm_mtz:
        logger.info('unmerged mtz exists')
#        if not os.path.isfile(unm_name):
#            logger.error('/bin/cp {0!s} .'.format(unm_mtz))
#            os.system('/bin/cp {0!s} .'.format(unm_mtz))
#            run_mrfana(logger, unm_name)
    # taking out this line; will keep entire cif file, not just header
#    if not os.path.isfile(cif_name):
    if not os.path.isfile('process_header.cif'):
        write_mmcif_header(logger, cif, cif_name, collection_date, wavelength)
    create_process_symlink(mtz_name, log_name, cif_name)
    mtz = os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run), pipeline,
                       'process.mtz')
    log = os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run), pipeline,
                       'process.log')
    cif = os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run), pipeline,
                       'process.cif')
    return mtz, log, cif, mrfana_ciffile

def check_first_line_of_aimless_mrfana_cif(logger, mrf_name, projectDir, sample, proposal, session, run, pipeline):
    os.chdir(os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run), pipeline))
    logger.info('checking start line of mrfana...')
    with open(mrf_name, 'r+') as f:
        content = f.read()
        if not content.startswith('data_mrfana'):
            logger.info('inserting data_mrfana as first line...')
            f.seek(0, 0)
            f.write('data_mrfana\n' + content)

def run_mrfana(logger, unm_name):
    if os.path.isfile('unmerged_total.cif'):
        logger.info('file exists: unmerged_total.cif; skipping...')
    else:
        logger.info('running MRFANA on UNMERGED MTZ')
        cmd = (
            'module load gopresto BUSTER\n'
            'mrfana {0!s} -cif unmerged.cif'.format(unm_name)
        )
        logger.info('command: {0!s}'.format(cmd))
        os.system(cmd)
        if os.path.isfile('unmerged_total.cif'):
            logger.info('ran mrfana successfully')
            check_first_line_of_mrfana_cif(logger)


def check_first_line_of_mrfana_cif(logger):
    logger.info('checking start line of mrfana...')
    with open('unmerged_total.cif', 'r+') as f:
        content = f.read()
        if not content.startswith('data_mrfana'):
            logger.info('inserting data_mrfana as first line...')
            f.seek(0, 0)
            f.write('data_mrfana\n' + content)


def create_process_symlink(mtz_name, log_name, cif_name):
    if not os.path.isfile('process.mtz'):
        os.system('ln -s {0!s} process.mtz'.format(mtz_name))
    if not os.path.isfile('process.log'):
        os.system('ln -s {0!s} process.log'.format(log_name))
    if not os.path.isfile('process.cif'):
        os.system('ln -s {0!s} process.cif'.format(cif_name))
#        os.system('ln -s {0!s} process.cif'.format(cif_name.replace('.bz2', '')))


def mtz_info(mtzfile):
    mtzDict = {}
    mtz = gemmi.read_mtz_file(mtzfile)
    mtzDict['unitcell'] = str(mtz.cell.parameters).replace('(', '').replace(')', '')
    mtzDict['a'] = mtz.cell.a
    mtzDict['b'] = mtz.cell.b
    mtzDict['c'] = mtz.cell.c
    mtzDict['alpha'] = mtz.cell.alpha
    mtzDict['beta'] = mtz.cell.beta
    mtzDict['gamma'] = mtz.cell.gamma
    mtzDict['unitcell_volume'] = mtz.cell.volume
    mtzDict['point_group'] = mtz.spacegroup.point_group_hm()
    mtzDict['resolution_high'] = mtz.resolution_high()
    mtzDict['space_group'] = mtz.spacegroup.hm
    mtzDict['wavelength'] = mtz.dataset(1).wavelength
#    mtzDict['lattice'] = mtz.spacegroup.crystal_system_str()
    mtzDict['lattice'] = mtz.spacegroup.hm[0]
    return mtzDict

def update_mrfana_cif(logger, ciffile):
    logger.warning('updating mrfana20 cif file...')
    os.chdir(ciffile[:ciffile.rfind('/')])
    out = "data_mrfana\n"
    found = False
    for line in open(os.path.realpath(ciffile)):
        if line.startswith('_reflns.details'):
            out += '_reflns.details ?\n'
            continue
        if line.startswith(';') and not found:
            found = True
            continue
        if line.startswith(';') and found:
            found = False
            continue
        if not found:
            out += line

    f = open('aimless.mrfana20.edited.cif', 'w')
    f.write(out)
    f.close()
    if os.path.isfile('process.cif'):
        os.system('/bin/rm process.cif')
    os.system('ln -s aimless.mrfana20.edited.cif process.cif')
    logger.warning('finished updating mrfana20 cif file')

def make_cif_header_only(logger, ciffile):
    doc = None
    os.chdir(ciffile[:ciffile.rfind('/')])
    logger.info('trying to read process_header.cif instead')
    if os.path.isfile('process_header.cif'):
        logger.info('reading process_header.cif')
        doc = gemmi.cif.read_file('process_header.cif')
    return doc

def cif_info(logger, ciffile):
    cifDict = {}
    try:
#        if ciffile:
#            proc_header = ciffile.replace('process.cif', 'process_header.cif')
#        else:
#            logger.error('ciffile is None; skipping...')
#            return cifDict
#        if os.path.isfile(proc_header):
#            logger.info(f"trying to read {proc_header} with gemmi")
#            doc = gemmi.cif.read_file(proc_header)
#        else:
        logger.info(f"trying to read {ciffile} with gemmi")
        doc = gemmi.cif.read_file(ciffile)
    except ValueError:
        logger.error('gemmi throws a ValueError for {0!s}'.format(ciffile))
        if "mrfana20" in os.path.realpath(ciffile):
            logger.warning('this looks like a file from mrfana...')
            update_mrfana_cif(logger, ciffile)
            doc = gemmi.cif.read_file(ciffile)
        else:
            logger.warning(f'there is something wrong with {ciffile}')
            newciffile = ciffile.replace('process.cif', 'aimless.mrfana20.cif')
#            newciffile = ciffile.replace('process.cif', 'HDF5_1//aimless.mrfana20.cif')
            logger.info(f'checking if mrfana cif file exisits: {newciffile}')
            if os.path.isfile(newciffile):
                logger.info('file exists; trying to read it...')
                doc = gemmi.cif.read_file(newciffile)
#            doc = make_cif_header_only(logger, ciffile)
            if not doc:
                logger.error('could not read mrfana cif either')
#                logger.error('could not read cif header either')
                return cifDict
    logger.info('iterating through doc...')
    for block in doc:
#        if block.find_pair('_symmetry.space_group_name_H-M'):
#            cifDict['space_group'] = str(block.find_pair('_symmetry.space_group_name_H-M')[1])
#        if block.find_pair('_cell.length_a'):
#            cifDict['a'] = str(block.find_pair('_cell.length_a')[1])
#            cifDict['b'] = str(block.find_pair('_cell.length_b')[1])
#            cifDict['c'] = str(block.find_pair('_cell.length_c')[1])
#            cifDict['alpha'] = str(block.find_pair('_cell.angle_alpha')[1])
#            cifDict['beta'] = str(block.find_pair('_cell.angle_beta')[1])
#            cifDict['gamma'] = str(block.find_pair('_cell.angle_gamma')[1])
        if block.find_pair('_reflns.d_resolution_low'):
            cifDict['reso_low'] = str(round(float(block.find_pair('_reflns.d_resolution_low')[1]), 2))
#            cifDict['reso_low'] = str(block.find_pair('_reflns.d_resolution_low')[1])
            cifDict['reso_high'] = str(round(float(block.find_pair('_reflns.d_resolution_high')[1]), 2))
#        if block.find_pair('_reflns.pdbx_netI_over_sigmaI'):
            cifDict['pdbx_netI_over_sigmaI'] = str(block.find_pair('_reflns.pdbx_netI_over_sigmaI')[1])
            cifDict['pdbx_redundancy'] = str(block.find_pair('_reflns.pdbx_redundancy')[1])
            if 'staraniso' in ciffile:
                try:
                    cifDict['percent_possible_obs'] = str(block.find_pair('_reflns.pdbx_percent_possible_spherical')[1])
                except TypeError:
                    pass
            elif 'xia2' in ciffile:
                cifDict['percent_possible_obs'] = str(round((float(block.find_pair('_reflns.percent_possible_obs')[1])*100), 1))
            else:
                cifDict['percent_possible_obs'] = str(block.find_pair('_reflns.percent_possible_obs')[1])
            cifDict['pdbx_number_measured_all'] = str(block.find_pair('_reflns.pdbx_number_measured_all')[1])
            cifDict['pdbx_CC_half'] = str(block.find_pair('_reflns.pdbx_CC_half')[1])

        if block.find_loop('_reflns_shell.Rmerge_I_obs'):
            Rmerge_I_obs = list(block.find_loop('_reflns_shell.Rmerge_I_obs'))
            cifDict['Rmerge_I_obs_low'] = str(round(float(Rmerge_I_obs[0]), 2))

        if block.find_loop('_reflns_shell.meanI_over_sigI_obs'):
            meanI_over_sigI_obs = list(block.find_loop('_reflns_shell.meanI_over_sigI_obs'))
            cifDict['meanI_over_sigI_obs_high'] = str(round(float(meanI_over_sigI_obs[len(meanI_over_sigI_obs)-1]), 2))
    return cifDict


def get_status(logger, mtzfile, mtz, ciffile, status, mrfana_ciffile, proc_header):
#    if mtzfile and (ciffile or mrfana_ciffile):
    if mtzfile and os.path.isfile(proc_header):
        mtzDict = mtz_info(mtzfile)
        if float(mtzDict['resolution_high']) < 2.5:
            status = 'OK'
            if os.path.isfile(proc_header):
                cif = cif_info(logger, proc_header)
            elif mrfana_ciffile:
                cif = cif_info(logger, mrfana_ciffile)
            if 'Rmerge_I_obs_low' in cif:
                if float(cif['Rmerge_I_obs_low']) > 0.15:
                    logger.error('Rmerge of {0!s} is higher than 15%'.format(cif['Rmerge_I_obs_low']))
                    status = 'FAIL - high Rmerge (low)'
            else:
                status = 'FAIL - unknown'
        elif float(mtzDict['resolution_high']) >= 2.8 and float(mtzDict['resolution_high']) < 3.2:
            status = 'FAIL - medium resolution'
        else:
            status = 'FAIL - low resolution'
    else:
        status = 'FAIL - no processing result'
    if mtzfile:
        status = check_if_salt_lattice(logger, mtz, status)
    return status


def check_if_salt_lattice(logger, mtz, status):
    if float(mtz['a']) < 20 or float(mtz['b']) < 20 or float(mtz['c']) < 20:
        logger.warning('at least one unit cell axis is shorter than 20A (see below); looks like salt or compound...')
        logger.info('unit cell: {0!s}'.format(mtz['unitcell']))
        status = 'FAIL - SALT?!'
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
#    if pipeline:
#        os.chdir(pipeline)
    d = {}
    d['sample'] = sample
    d['collection_date'] = collection_date
    d['run'] = run
    d['proposal'] = proposal
    d['session'] = session
    d['protein'] = protein
    d['status'] = ""
    d['master'] = master
    d['pipeline'] = ""
    logger.info('saving info.json file without pipeline and status information...')
    with open("info.json", "w") as outfile:
        json.dump(d, outfile)
    if pipeline:
        logger.info('changing to pipeline folder')
        os.chdir(pipeline)
        d['pipeline'] = pipeline
        d['status'] = status
        logger.info('saving info.json file (again)...')
        with open("info.json", "w") as outfile:
            json.dump(d, outfile)

def read_info_json_file(logger, json_file):
    proposal = None
    session = None
    run = None
    logger.info("reading {0!s}".format(json_file))
    d = None
    with open(json_file) as jfile:
        d = json.load(jfile)
    if d:
        proposal = d['proposal']
        session = d['session']
        run = d['run']
        collection_date = d['collection_date']
    return proposal, session, run, collection_date

def read_reference_pdb_files(logger, projectDir):
    logger.info('reading point group, lattice and unit cell volume from reference pdb files...')
    os.chdir(os.path.join(projectDir, '0-model'))
    ref_dict = {}
    for pdbfile in glob.glob('*.pdb'):
        logger.info('current pdbfile: {0!s}'.format(pdbfile))
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


#def retain_results_with_similar_ucvol_and_pg_as_ref_pdb(logger, proc_dict, ref_dict):
#    logger.info('checking if PDB files in 0-model folder are similar to the processing MTZ files...')
#    match_dict = {}
#    if ref_dict:
#        for f in proc_dict:
#            logger.info('current CIF ' + f)
#            pgr_mtz = proc_dict[f]['point_group']
#            ucv_mtz = float(proc_dict[f]['unitcell_volume'])
#            lat_mtz = proc_dict[f]['lattice']
#            for p in ref_dict:
#                pgr_pdb = ref_dict[p][0]
#                ucv_pdb = float(ref_dict[p][2])
#                lat_pdb = ref_dict[p][1]
#                ucv_diff = abs((ucv_mtz-ucv_pdb))/ucv_pdb
#                if pgr_mtz == pgr_pdb and lat_mtz == lat_pdb and ucv_diff < 0.1:
#                    logger.info('lattice, point group and unit cell volume of MTZ file matches {0!s}'.format(p))
#                    match_dict[f] = proc_dict[f]
#                    break
#            logger.warning('no match found!')
#    else:
#        logger.warning('seems like no reference PDB files were provided; skipping this selection step...')
#    if match_dict:
#        proc_dict = match_dict
#    else:
#        logger.warning('none of the PDB files appears to be similar to any of the auto-processing MTZ files')
#    return proc_dict

def retain_results_with_similar_ucvol_and_pg_as_ref_pdb(logger, proc_list, ref_dict):
    logger.info('checking if PDB files in 0-model folder are similar to the processing MTZ files...')
    match_list = []
    match_found = False
    if ref_dict:
        for d in proc_list:
            logger.info('current pipeline {0!s}'.format(d['autoproc_pipeline']))
#            pgr_mtz = proc_dict[f]['point_group']
            pgr_mtz = d['sym_point_group']
#            ucv_mtz = float(proc_dict[f]['unitcell_volume'])
            try:
                ucv_mtz = float(d['cell_volume'])
            except ValueError:
                logger.error('cannot convert unit cell volume to float: {0!s}'.format(d['cell_volume']))
                continue
#            lat_mtz = proc_dict[f]['lattice']
            lat_mtz = d['sym_lattice']
            logger.info('CIF -> lat: {0!s} - pg: {1!s} - ucv: {2!s}'.format(lat_mtz, pgr_mtz, ucv_mtz))
            for p in ref_dict:
                pgr_pdb = ref_dict[p][0]
                ucv_pdb = float(ref_dict[p][2])
                lat_pdb = ref_dict[p][1]
                logger.info('CIF -> lat: {0!s} - pg: {1!s} - ucv: {2!s}'.format(lat_pdb, pgr_pdb, ucv_pdb))
                ucv_diff = abs((ucv_mtz-ucv_pdb))/ucv_pdb
                logger.info('unit cell volume difference = {0!s}'.format(ucv_diff))
                logger.info('unit cell volume difference = {0!s} %'.format(round(ucv_diff*100,2)))
                if pgr_mtz == pgr_pdb and lat_mtz == lat_pdb and ucv_diff < 0.1:
                    logger.info('lattice, point group and unit cell volume of MTZ file matches {0!s}'.format(p))
                    match_list.append(d)
#                    match_dict[f] = proc_dict[f]
                    break
                else:
                    logger.warning('no match found!')
    else:
        logger.warning('seems like no reference PDB files were provided; skipping this selection step...')
#    if match_dict:
#        proc_dict = match_dict
    if match_list:
        proc_list = match_list
        match_found = True
    else:
        logger.warning('none of the PDB files appears to be similar to any of the auto-processing MTZ files')
#        proc_list = []
        logger.info('will carry over all valid auto-processing files to the next stage')
    return proc_list, match_found


#def retain_results_with_good_low_reso_rmerge(logger, proc_dict):
#    logger.info('checking of auto-processing MTZ files have low resolution Rmerge values below 10%...')
#    match_dict = {}
#    for f in proc_dict:
#        if float(proc_dict[f]['Rmerge_I_obs_low']) < max_allowed_Rmerge_I_obs_low():
#            logger.info('{0!s} - Rmerge(low): {1!s}'.format(f, proc_dict[f]['Rmerge_I_obs_low']))
#            match_dict[f] = proc_dict[f]
#        else:
#            logger.error('{0!s} - Rmerge(low): {1!s}'.format(f, proc_dict[f]['Rmerge_I_obs_low']))
#    if match_dict:
#        proc_dict = match_dict
#    else:
#        logger.error('did not find any MTZ file with Rmerge (low) below {0!s}; skipping sample...'.format(
#            max_allowed_Rmerge_I_obs_low()))
#        proc_dict = {}
#    return proc_dict


def retain_results_with_good_low_reso_rmerge(logger, proc_list):
    logger.info('checking of auto-processing MTZ files have low resolution Rmerge values below 10%...')
    match_list = []
    for d in proc_list:
        logger.info('current pipeline: {0!s}'.format(d['autoproc_pipeline']))
        try:
            if float(d['reflns_inner_pdbx_Rmerge_I_obs']) < max_allowed_Rmerge_I_obs_low():
#            logger.info('{0!s} - Rmerge(low): {1!s}'.format(f, proc_dict[f]['Rmerge_I_obs_low']))
                logger.info('acceptable low resolution Rmerge: {0!s}'.format(d['reflns_inner_pdbx_Rmerge_I_obs']))
                match_list.append(d)
            else:
#            logger.error('{0!s} - Rmerge(low): {1!s}'.format(f, proc_dict[f]['Rmerge_I_obs_low']))
                logger.warning('low resolution Rmerge is too high: {0!s}'.format(d['reflns_inner_pdbx_Rmerge_I_obs']))
        except ValueError:
            logger.error(f"value for reflns_inner_pdbx_Rmerge_I_obs cannot be converted to string: {d}")
            logger.warning("removing item from proc_list")
            proc_list.remove(d)
    if match_list:
        proc_list = match_list
    else:
#        logger.error('did not find any MTZ file with Rmerge (low) below {0!s}; skipping sample...'.format(
#            max_allowed_Rmerge_I_obs_low()))
#        proc_list = []
        logger.warning('did not find any MTZ file with Rmerge (low) below {0!s}; will take it anyway...'.format(
            max_allowed_Rmerge_I_obs_low()))
        logger.info('will still cary over all valid files...')
        logger.info(proc_list)
    return proc_list


#def retain_results_which_fit_selection_criterion(logger, proc_dict, select_criterion):
#    logger.info('selecting auto-processing results based on {0!s}...'.format(select_criterion))
#    match_list = []
#    backup_list = []
#    found_selected_pipeline = False
#    for f in proc_dict:
#        reso_high = proc_dict[f]['reso_high']
#        if select_criterion.startswith('reso'):
#            logger.info('added {0!s} with high resolution limit of {1!s} A'.format(f, reso_high))
#            match_list.append([f, float(reso_high)])
#        elif proc_dict[f]['pipeline'] == select_criterion:
#            logger.info('added {0!s} with high resolution limit of {1!s} A'.format(f, reso_high))
#            found_selected_pipeline = True
#            match_list.append([f, float(reso_high)])
#        else:
#            logger.warning('MTZ does not match criteria, but added {0!s} with high resolution limit of {1!s} A'.format(
#                f, reso_high))
#            backup_list.append([f, float(reso_high)])
#    if not match_list:
#        logger.warning('none of the MTZ files fulfilled the selection criteria, but will select the one with highest resolution')
#        match_list = backup_list
#    if proc_dict:
#        logger.info('current list of matching auto-processing results: {0!s}'.format(match_list))
#        bestcif = min(match_list, key=lambda x: x[1])[0]
#        logger.info('--> {0!s}'.format(bestcif))
#    else:
#        bestcif = None
#    return bestcif, found_selected_pipeline

def retain_results_which_fit_selection_criterion(logger, dal, proc_list, select_criterion):
    logger.info('selecting auto-processing results based on: {0!s}'.format(select_criterion))
    match_list = []
    backup_list = []
    found_selected_pipeline = False
    for d in proc_list:
        logger.info(f"current pipeline: {d['autoproc_pipeline']}")
        processing_id = d['processing_id']
        processing_outcome = "unknown"
        reso_high = d['reflns_d_resolution_high']
        if select_criterion.startswith('reso'):
            logger.info('added {0!s} with high resolution limit of {1!s} A'.format(d['autoproc_pipeline'], reso_high))
            processing_outcome = "success - is selected highres"
            match_list.append([d, float(reso_high)])
        elif select_criterion.lower() == 'autoproc' and d['autoproc_pipeline'].lower() == 'autoproc' and d['staraniso'] == False:
            logger.info('added {0!s} with high resolution limit of {1!s} A'.format(d['autoproc_pipeline'], reso_high))
            found_selected_pipeline = True
            processing_outcome = "success - fits selected pipeline"
            match_list.append([d, float(reso_high)])
        elif select_criterion.lower() == 'staraniso' and d['autoproc_pipeline'].lower() == 'autoproc' and d['staraniso'] == True:
            logger.info('added {0!s} with high resolution limit of {1!s} A'.format(d['autoproc_pipeline'], reso_high))
            found_selected_pipeline = True
            processing_outcome = "success - fits selected pipeline"
            match_list.append([d, float(reso_high)])
        else:
            logger.warning('MTZ does not match criteria, but added {0!s} with high resolution limit of {1!s} A'.format(
                d['autoproc_pipeline'], reso_high))
            processing_outcome = "fail - not selected pipeline"
            try:
                backup_list.append([d, float(reso_high)])
            except ValueError:
                logger.error('seems like that the high resolution limit in the database cannot be converted to float')
        processdb.update_processing_outcome(logger, dal, processing_id, processing_outcome)
    if not match_list:
        logger.warning('none of the MTZ files fulfilled the selection criteria, but will select the one with highest resolution')
        match_list = backup_list
    if proc_list:
#        logger.info('current list of matching auto-processing results: {0!s}'.format(match_list))
        best = min(match_list, key=lambda x: x[1])[0]
        logger.info('best match --> {0!s}'.format(best['autoproc_pipeline']))
    else:
        best = None
    logger.warning('1 -> {0!s}'.format(best))
    return best, found_selected_pipeline

#def get_result_with_highest_resolution(logger, dal, proc_list, select_criterion):
#    match_list = []
#    for d in proc_list:
#        processing_id = d['processing_id']
#        processing_outcome = "unknown"
#        reso_high = d['reflns_d_resolution_high']
#        match_list.append([d, float(reso_high)])
#
#        if select_criterion.startswith('reso'):
#            logger.info('added {0!s} with high resolution limit of {1!s} A'.format(d['autoproc_pipeline'], reso_high))
#            processing_outcome = "success - is selected highres"


def check_if_best_result_is_from_select_pipeline(logger, sample, found_selected_pipeline, not_fitting_pipeline_list, select_criterion):
#    pipeline_list = ['xia2dials', 'autoproc', 'xia2xds', 'staraniso']
    pipeline_list = ['xia2', 'autoproc', 'staraniso']
    if select_criterion in pipeline_list:
        if not found_selected_pipeline:
            logger.warning('found a reasonable data processing result, but not for the selected auto-processing pipeline')
            not_fitting_pipeline_list.append(sample)
    return not_fitting_pipeline_list

def report_not_fitting_pipelines(logger, not_fitting_pipeline_list, processDir, projectDir, fragmaxcsv):
    csv_out = ''
    logger.warning('the following samples do not have a matching result from the selected auto-processing pipeline:')
    for sample in not_fitting_pipeline_list:
        logger.info(' --> {0!s}'.format(sample))
        csv_out += sample + ',\n'
    save_reprocess_csv_file(logger, csv_out, processDir, projectDir, fragmaxcsv)

def save_reprocess_csv_file(logger, csv_out, processDir, projectDir, fragmaxcsv):
    logger.info('saving reprocess.csv file as ')
    os.chdir(os.path.join(projectDir, 'script'))
    f = open('reprocess.csv', 'w')
    f.write(csv_out)
    f.close()
    suggest_reprocessing_input(logger, processDir, projectDir, fragmaxcsv, "reprocess.csv")

def suggest_reprocessing_input(logger, processDir, projectDir, fragmaxcsv, reprocesscsv):
    cmd = (
        'python /data/staff/biomax/tobias/software/MAXIV_tools/1-process.py '
        '-i {0!s} '.format(processDir) +
        '-o {0!s} '.format(projectDir) +
        '-f {0!s} '.format(fragmaxcsv) +
        '-r {0!s}'.format(reprocesscsv)
    )
    logger.info('try running the following command to reprocess the datasets with the selected pipeline:')
    logger.info(cmd)

#def link_process_results(logger, projectDir, sample, bestcif):
#    logger.info('creating symlinks in {0!s}'.format(os.path.join(projectDir, "1-process", sample)))
#    os.chdir(os.path.join(projectDir, "1-process", sample))
#    pipeline = bestcif.split('/')[1]
##    json_info = bestcif.replace(pipeline +'/process.cif', 'info.json')
#    json_info = bestcif.replace('process.cif', 'info.json')
#    if not os.path.isdir('process.mtz'):
#        os.system('ln -s {0!s} .'.format(bestcif.replace('.cif', '.mtz')))
#    if not os.path.isdir('process.log'):
#        os.system('ln -s {0!s} .'.format(bestcif.replace('.cif', '.log')))
#    if not os.path.isdir('process.cif'):
#        os.system('ln -s {0!s} .'.format(bestcif))
#    if not os.path.isdir('info.json'):
#        os.system('ln -s {0!s} .'.format(json_info))


def link_process_results(logger, projectDir, sample, best):
    logger.info('creating symlinks in {0!s}'.format(os.path.join(projectDir, "1-process", sample)))
    logger.warning('2 -> {0!s}'.format(best))
    os.chdir(os.path.join(projectDir, "1-process", sample))
#    pipeline = bestcif.split('/')[1]
#    json_info = bestcif.replace(pipeline +'/process.cif', 'info.json')
#    json_info = bestcif.replace('process.cif', 'info.json')
    json_info = os.path.relpath(best['processing_cif_file']).replace('process.cif', 'info.json')
    logger.info('>>> MTZ {0!s}'.format(best['processing_mtz_file']))
    logger.info('>>> LOG {0!s}'.format(best['processing_log_file']))
    logger.info('>>> CIF {0!s}'.format(best['processing_cif_file']))
    if not os.path.isdir('process.mtz'):
#        os.system('ln -s {0!s} .'.format(bestcif.replace('.cif', '.mtz')))
    # note: need to wrap path into realpath/ relpath, because,
    # /data/visitors/... really is /gpfs/offline/visitors and the relpath statement alone will turn
    # into ../../../../../../../../data/... but once a realpath come first it will be ./...
        logger.info('ln -s {0!s} process.mtz'.format(os.path.relpath(os.path.realpath(best['processing_mtz_file']))))
        os.system('ln -s {0!s} process.mtz'.format(os.path.relpath(os.path.realpath(best['processing_mtz_file']))))
    if not os.path.isdir('process.log'):
#        os.system('ln -s {0!s} .'.format(bestcif.replace('.cif', '.log')))
        os.system('ln -s {0!s} process.log'.format(os.path.relpath(os.path.realpath(best['processing_log_file']))))
    if not os.path.isdir('process.cif'):
#        os.system('ln -s {0!s} .'.format(bestcif))
        os.system('ln -s {0!s} process.cif'.format(os.path.relpath(os.path.realpath(best['processing_cif_file']))))
    if not os.path.isdir('info.json'):
        os.system('ln -s {0!s} .'.format(os.path.relpath(os.path.realpath(json_info))))

def link_info_json_file(logger, projectDir, sample):
    foundDozor = False
    os.chdir(os.path.join(projectDir, "1-process", sample))
    for f in glob.glob(os.path.join('*', 'info.json')):
        logger.info('linking {0!s} to sample directory'.format(f))
        os.system('ln -s {0!s} .'.format(f))
        foundDozor = True
    if not foundDozor:
        # this situation may occur if data processing timed out
        for f in glob.glob(os.path.join('*', '*', 'info.json')):
            logger.info('linking {0!s} to sample directory'.format(f))
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


#def read_data_collection_stats(logger, ciffile, proc_dict):
#    subfolder = ciffile.split('/')[0]
#    pipeline = ciffile.split('/')[1]
#    mtzfile = os.path.join(ciffile.replace('.cif', '.mtz'))
#    if os.path.isfile(mtzfile):
#        logger.info('found {0!s} file'.format(mtzfile))
#        mtz = mtz_info(mtzfile)
#        mtz['pipeline'] = pipeline
#        mtz['subfolder'] = subfolder
#        cif = cif_info(logger, ciffile)
#        mtz.update(cif)
#        proc_dict[os.path.join(ciffile)] = mtz
#    else:
#        logger.error('process.mtz does not exist in folder')
#    return proc_dict


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


def start_reprocessing(logger):
    logger.info('===================================================================================')
    logger.info('>>> START: reprocessing selected files')
    logger.info('===================================================================================')


def end_reprocessing(logger):
    logger.info('===================================================================================')
    logger.info('>>> END: finished reprocessing selected files')
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
#    supported_options = ['xia2dials', 'xia2xds', 'autoproc', 'staraniso', 'resolution']
    supported_options = ['xia2', 'autoproc', 'staraniso', 'resolution']
    if select_criterion in supported_options:
        logger.info('option exists')
    else:
        logger.error('ERROR: option does not exist')
        passed = False
    return passed

def check_if_database_file_exists(logger, db_file, passed):
    logger.info('checking db_file: {0!s}'.format(db_file))
    if os.path.isfile(db_file):
        logger.info('db_file exists')
    else:
        logger.error('db_file does not exist')
        passed = False
    return passed


def run_checks(logger, processDir, projectDir, fragmaxcsv, select, select_criterion, db_file):
    logger.info('checking input file and command line options...')
    passed = check_if_process_directory_exists(logger, processDir, select)
    passed = check_if_project_directory_exists(logger, projectDir, passed)
    passed = check_if_fragmaxcsv_exists(logger, fragmaxcsv, passed)
    passed = check_pdbfiles_in_pdbdir(logger, projectDir, passed)
    passed = check_select_criterion(logger, select_criterion, passed)
    passed = check_if_database_file_exists(logger, db_file, passed)
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
        q = input("\n>>> Do you want to continue? (y/n) ")
    if not q.lower() == 'y':
        logger.info('you chose not to continue at this point; exciting program...')
        sys.exit(2)


def get_proc_dict():
    proc_dict = {
        'pipeline': '',
        'space_group': '',
        'unit_cell': ''
    }
    return proc_dict


def ask_for_spg_and_unit_cell(logger):
    proc_dict = get_proc_dict()
    q = ''
    if sys.version[0] == '2':
        q = raw_input("\n>>> specify processing pipeline: ")
    else:
        q = input("\n>>> specify processing pipeline: ")
    if q.replace(' ','') == '':
        logger.info('you must choose a processing pipeline; exciting program...')
        sys.exit(2)
    proc_dict['pipeline'] = q.replace(' ','')
    q = ''
    if sys.version[0] == '2':
        q = raw_input("\n>>> specify space group or press return: ")
    else:
        q = input("\n>>> specify space group or press return: ")
    proc_dict['space_group'] = q.replace(' ','')
    q = ''
    if sys.version[0] == '2':
        q = raw_input("\n>>> specify unit cell or press return: ")
    else:
        q = input("\n>>> specify unit cell or press return: ")
#    proc_dict['unit_cell'] = q.replace(' ','')
    proc_dict['unit_cell'] = ','.join(q.split())
    return proc_dict


def get_script_dict(pipeline, n_jobs, now):
    cmd = maxiv_header(pipeline)
    cmd += modules_to_load(pipeline)
    script_dict = {}
    for i in range(n_jobs):
        script_dict[pipeline+'_{0!s}_{1!s}.sh'.format(now, i)] = cmd
#    script_dict[pipeline+'_1.sh'] = cmd
    return script_dict


def modules_to_load(pipeline):
    module = ''
    if pipeline.startswith('xia2'):
        module = 'module load gopresto CCP4\n'
    elif pipeline.startswith('autoproc'):
#        module = 'module load gopresto BUSTER\n'
        module = 'module load gopresto autoPROC/20240123-PReSTO-10.0\n'
#        module = 'module load gopresto XDS\n'
#        module += 'module load gopresto CCP4\n'
#        module += 'source /mxn/groups/sw/mxsw/autoPROC_20240123/setup.sh\n'
    return module


def maxiv_header(pipeline):
    header = (
        '#!/bin/bash\n'
        '#SBATCH --time=99:00:00\n'
        '#SBATCH --job-name={0!s}\n'.format(pipeline) +
        '#SBATCH --cpus-per-task=48\n'
        '#SBATCH --exclusive\n'
        '#SBATCH --nodes=1\n'
    )
    return header


def pipeline_cmd(pipeline, proc_dict, master_file):
    extra_cmd_xia = ''
    extra_cmd_autoproc = ''
    if proc_dict['space_group']:
        extra_cmd_xia += 'space_group=' + proc_dict['space_group'] + ' '
        extra_cmd_autoproc +=  'symm="{0!s}"'.format(proc_dict['space_group']) + ' '
    if proc_dict['unit_cell'] and proc_dict['unit_cell'] != "":
        extra_cmd_xia += 'unit_cell=' + proc_dict['unit_cell'] + ' '
        extra_cmd_autoproc += 'cell="{0!s}"'.format(proc_dict['unit_cell'].replace(',', ' ')) + ' '
#    print('pipeline', pipeline)
    cmd = ''
    if pipeline.startswith('xia2dials'):
        cmd = 'xia2 pipeline=dials image={0!s} {1!s}'.format(master_file, extra_cmd_xia)
    elif pipeline.startswith('xia2xds'):
        cmd = 'xia2 pipeline=3dii image={0!s} {1!s}'.format(master_file, extra_cmd_xia)
    elif pipeline.startswith('autoproc'):
        cmd = 'process -h5 {0!s} -nthreads 48 {1!s} '.format(master_file, extra_cmd_autoproc)
    return cmd


def get_proc_folder(projectDir, sample, proposal, session, run, pipeline):
    p = os.path.join(projectDir, '1-process', sample, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run), pipeline)
    return p


def add_cmd_to_script_dict(logger, script_dict, counter, pipeline, proc_dict, proc_folder, master_file, now):
    logger.info('saving shell scripts for manual auto-processing...')
    script_dict[pipeline + '_{0!s}_{1!s}.sh'.format(now, counter)] += 'cd ' + proc_folder + '\n'
    script_dict[pipeline + '_{0!s}_{1!s}.sh'.format(now, counter)] += pipeline_cmd(pipeline, proc_dict, master_file) + '\n'
    return script_dict


def save_proc_scripts(logger, projectDir, script_dict):
    os.chdir(os.path.join(projectDir, 'tmp'))
    for script in script_dict:
        f = open(script, 'w')
        f.write(script_dict[script])
        f.close()

def check_if_to_annotate_and_reprocess(logger, missing_dict, dal):
    if sys.version[0] == '2':
        q = raw_input("\n>>> Do you want to annotate the affected samples? (y/n) ")
    else:
        q = input("\n>>> Do you want to annotate the affected samples? (y/n) ")
    if not q.lower() == 'y':
        logger.info('you chose not to manually annotate the affected samples!')
        logger.info('OK, then we are done here for the moment... bye, bye!')
    else:
        annotate_failed_datasets(logger, missing_dict, dal)

def get_fail_dict(logger):
    fail_dict = {
        '1': 'fail - no diffraction',
        '2': 'fail - weak diffraction',
        '3': 'fail - misaligned',
        '4': 'fail - loop broke',
        '5': 'fail - loop empty',
        '6': 'fail - no X-rays',
        '7': 'fail - salt/ compound crystal',
        '8': 'fail - data processing',
        '9': 'fail - no matching model',
        '10': 'fail - tiny crystal',
        '11': 'fail - queue stopped',
        '12': 'fail - unknown'
    }
    show_scoring_option(logger, fail_dict)
    return fail_dict

def get_accepted_scores(logger):
    fail_dict = get_fail_dict(logger)
    accepted_scores = []
    for n in fail_dict:
        accepted_scores.append(n)
    return accepted_scores

def show_scoring_option(logger, fail_dict):
    logger.info('scores - overview')
    for i in fail_dict:
        logger.info(' -> {0!s} == {1!s}'.format(i, fail_dict[i]))

def enter_score(logger, dal, dataset, fail_dict):
    sample = dataset[0]
    proposal = dataset[1]
    session = dataset[2]
    accepted_scores = get_accepted_scores(logger)
    if sys.version[0] == '2':
        q = raw_input(">>> {0!s}: ".format(sample))
    else:
        q = input(">>> {0!s}: ".format(sample))
    if not q in accepted_scores:
        logger.error('entry not in allowed scores; skipping...')
    else:
        mounted_crystal_id = processdb.get_mounted_crystal_id(dal, sample)
        d = {
            'data_collection_outcome': fail_dict[q]
        }
        w = {
            'mounted_crystal_id': mounted_crystal_id,
            'mounted_crystal_code': sample,
            'proposal': proposal,
            'session': session,
            'is_dataset': True
        }
        processdb.update_xray_dataset_table_with_dataset_outcome(logger, dal, d, w)

def annotate_failed_datasets(logger, missing_dict, dal):
    fail_dict = get_fail_dict(logger)
    logger.info('here are samples where no datasets were collectied')
    for category in missing_dict:
        for dataset in missing_dict[category]:
            enter_score(logger, dal, dataset, fail_dict)

def select_to_reprocess_missing_datasets(logger, dataset, csv_out):
    sample = dataset[0]
    if sys.version[0] == '2':
        q = raw_input(">>> {0!s} (y/n): ".format(sample))
    else:
        q = input(">>> {0!s} (y/n): ".format(sample))
    if q.lower() == 'y':
        logger.info('you chose not to manually annotate the affected samples!')
        csv_out += sample + ',\n'
    return csv_out

def reprocess_missing_datasets(logger, missing_dict, processDir, projectDir, fragmaxcsv):
    logger.info('would you like to re-process some of the missing datasets?')
    csv_out = ''
    for dataset in missing_dict['mtz_file']:
        csv_out = select_to_reprocess_missing_datasets(logger, dataset, csv_out)
    if csv_out:
        logger.info('saving reprocess.csv file...')
        save_reprocess_csv_file(logger, csv_out, processDir, projectDir, fragmaxcsv)

