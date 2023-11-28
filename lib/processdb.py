import os.path

import gemmi

import sqlalchemy
from sqlalchemy.sql import select, func
from sqlalchemy import and_

from sqlalchemy.dialects import sqlite
#print(statement.compile(dialect=sqlite.dialect()))
import sys
import h5py

def get_mounted_crystal_id(dal, sample):
    q = select([dal.mounted_crystals_table.c.mounted_crystal_id]).where(
        dal.mounted_crystals_table.c.mounted_crystal_code == sample)
    rp = dal.connection.execute(q)
    result = rp.fetchall()
    return result[0][0]

def get_d_xray_dataset_table_dict(logger, dal, sample, proposal, session, beamline, run, create_date, master, dozor_plot, crystal_snapshot_list, foundDataset):
    logger.info('getting d_xray_dataset_table_dict for {0!s}'.format(sample))
    mounted_crystal_id = get_mounted_crystal_id(dal, sample)

    # read more information from hdf5 master file
    # https://docs.h5py.org/en/stable/quick.html

    d_xray_dataset_table_dict = {
        'mounted_crystal_id':   mounted_crystal_id,
        'mounted_crystal_code': sample,
        'beamline':             beamline,
        'proposal':             proposal,
        'session':              session,
        'run':                  run,
        'data_collection_date': create_date,
        'h5_master_file':   master
    }

    if dozor_plot:
        d_xray_dataset_table_dict['dozor_plot'] = dozor_plot

    key = ['crystal_snapshot_1', 'crystal_snapshot_2', 'crystal_snapshot_3', 'crystal_snapshot_3']
    if crystal_snapshot_list:
        for n, img in enumerate(crystal_snapshot_list):
            if n >= 4:
                logger.warning('too many crystal snapshots: {0!s} - {1!s}'.format(run, img))
            else:
                d_xray_dataset_table_dict[key[n]] = img

    d_xray_dataset_table_dict, foundDataset = read_master_file(logger, master, d_xray_dataset_table_dict, foundDataset)

    return d_xray_dataset_table_dict, foundDataset

def read_master_file(logger, master_file, d_xray_dataset_table_dict, foundDataset):
    logger.info('reading master .h5 file: {0!s}'.format(master_file))
    if master_file:
        f = h5py.File(master_file, 'r')
        #    list(f.keys()) there is most likely only 1 key 'entry'
        dset = f['entry']
    #    d_xray_dataset_table_dict['detector_distance'] = dset['instrument']['detector']['distance'].value
        d_xray_dataset_table_dict['detector_distance'] = dset['instrument']['detector']['distance'][()]
    #    d_xray_dataset_table_dict['wavelength'] = dset['sample']['beam']['incident_wavelength'].value
        d_xray_dataset_table_dict['wavelength'] = dset['sample']['beam']['incident_wavelength'][()]
        d_xray_dataset_table_dict['n_images'] = dset['sample']['goniometer']['omega'].shape[0]
        d_xray_dataset_table_dict['omega_range_total'] = dset['sample']['goniometer']['omega_range_total'][()]
        if d_xray_dataset_table_dict['omega_range_total'] > 30.0:
            d_xray_dataset_table_dict['is_dataset'] = True
            foundDataset = True
            d_xray_dataset_table_dict['data_collection_type'] = "rotation dataset"
            d_xray_dataset_table_dict['data_collection_outcome'] = "collected - processing pending"
        else:
            d_xray_dataset_table_dict['data_collection_type'] = "X-ray centring"
            d_xray_dataset_table_dict['data_collection_outcome'] = "X-ray centring"
    else:
        logger.error('master file does not exist!')
    return d_xray_dataset_table_dict, foundDataset

def insert_into_xray_dataset_table(logger, dal, d):
    logger.info('trying to insert into xray_dataset_table')
    try:
#        print(d)
        ins = dal.xray_dataset_table.insert().values(d)
#        print(ins.compile(dialect=sqlite.dialect()))
        dal.connection.execute(ins)
    except sqlalchemy.exc.IntegrityError as e:
        if "UNIQUE constraint failed" in str(e):
            logger.warning('entry exists (time soaked {0!s}); skipping'.format(d['mounted_crystal_code']))
        else:
            logger.error(str(e))

def update_xray_dataset_table_with_dataset_outcome(logger, dal, d, w):
    logger.info('updating xray_dataset table')
    u = dal.xray_dataset_table.update().values(d).where(and_(
        dal.xray_dataset_table.c.mounted_crystal_code == w['mounted_crystal_code'],
        dal.xray_dataset_table.c.mounted_crystal_id == w['mounted_crystal_id'],
        dal.xray_dataset_table.c.proposal == w['proposal'],
        dal.xray_dataset_table.c.session == w['session'],
        dal.xray_dataset_table.c.is_dataset == w['is_dataset']))
    dal.connection.execute(u)


def create_dummy_dataset_entry(logger, dal, sample, proposal, session):
    mounted_crystal_id = get_mounted_crystal_id(dal, sample)
    d = {
        'mounted_crystal_id': mounted_crystal_id,
        'mounted_crystal_code': sample,
        'proposal': proposal,
        'session': session,
        'run': 'dummy',
        'is_dataset': True,
        'data_collection_type': "dummy",
        'data_collection_outcome': "fail - no diffraction",
        'data_collection_comment': 'dummy entry - no dataset collected after x-ray alignment'
    }
    insert_into_xray_dataset_table(logger, dal, d)


def get_cell_sym_info(mtz, d):
    d['cell_length_a'] = mtz.cell.a
    d['cell_length_b'] = mtz.cell.b
    d['cell_length_c'] = mtz.cell.c
    d['cell_angle_alpha'] = mtz.cell.alpha
    d['cell_angle_beta'] = mtz.cell.beta
    d['cell_angle_gamma'] = mtz.cell.gamma
    d['cell_volume'] = mtz.cell.volume
    d['sym_lattice'] = mtz.spacegroup.hm[0]
    d['sym_point_group'] = mtz.spacegroup.point_group_hm()
    d['sym_lattice_point_group'] = str(d['sym_lattice']) + str(d['sym_point_group'])
    d['sym_space_group'] = mtz.spacegroup.hm
    d['sym_Int_Tables_number'] = mtz.spacegroup.number
    return d

def get_software_info(logger, block, d):
    data_reduction_software_list = ['XDS', 'DIALS']
    data_scaling_software_list = ['AIMLESS', 'DIALS']
    autoproc_pipeline_list = ['autoPROC', 'xia2']
    if block.find_loop('_software.name'):
        software = list(block.find_loop('_software.name'))
        version = list(block.find_loop('_software.version'))
        logger.info('software list: {0!s}'.format(list(block.find_loop('_software.name'))))
        for n,item in enumerate(software):
            if item in data_reduction_software_list and item != "STARANISO":
                d['data_reduction_software'] = item
                d['data_reduction_software_version'] = version[n]
            elif item == "STARANISO":
                d['staraniso'] = True
                d['staraniso_version'] = version[n]
            elif item in data_scaling_software_list:
                d['data_scaling_software'] = item
                d['data_scaling_software_version'] = version[n]
            elif item in autoproc_pipeline_list:
                d['autoproc_pipeline'] = item
                d['autoproc_pipeline_version'] = version[n]
    return d

def get_overall_stats(block, d):
    if block.find_pair('_reflns.d_resolution_low'):
        d['reflns_d_resolution_low'] = round(float(block.find_pair('_reflns.d_resolution_low')[1]), 2)
        d['reflns_d_resolution_high'] = round(float(block.find_pair('_reflns.d_resolution_high')[1]), 2)
        d['reflns_pdbx_netI_over_sigmaI'] = block.find_pair('_reflns.pdbx_netI_over_sigmaI')[1]
        d['reflns_pdbx_redundancy'] = block.find_pair('_reflns.pdbx_redundancy')[1]
        d['reflns_number_obs'] = block.find_pair('_reflns.number_obs')[1]
        d['reflns_percent_possible_obs'] = block.find_pair('_reflns.percent_possible_obs')[1]
        d['reflns_pdbx_Rmerge_I_obs'] = block.find_pair('_reflns.pdbx_Rmerge_I_obs')[1]
        d['reflns_pdbx_pdbx_Rrim_I_all'] = block.find_pair('_reflns.pdbx_Rrim_I_all')[1]
        d['reflns_pdbx_CC_half'] = block.find_pair('_reflns.pdbx_CC_half')[1]
    return d

def get_dataset_id(dal, mounted_crystal_code, proposal, session, run):
    q = select([dal.xray_dataset_table.c.dataset_id]).where(and_(
                dal.xray_dataset_table.c.mounted_crystal_code == mounted_crystal_code,
                dal.xray_dataset_table.c.proposal == proposal,
                dal.xray_dataset_table.c.session == session,
                dal.xray_dataset_table.c.run == run))
    rp = dal.connection.execute(q)
    r = rp.fetchall()
    idx = r[0][0]
    return idx

def get_lowres_stats(block, d):
    if block.find_loop('_reflns_shell.pdbx_ordinal'):
        if block.find_loop('_reflns_shell.d_res_high'):
            d['reflns_inner_d_resolution_high'] = list(block.find_loop('_reflns_shell.d_res_high'))[0]
        if block.find_loop('_reflns_shell.d_res_low'):
            d['reflns_inner_d_resolution_low'] = list(block.find_loop('_reflns_shell.d_res_low'))[0]
        if block.find_loop('_reflns_shell.number_measured_obs'):
            d['reflns_inner_number_obs'] = list(block.find_loop('_reflns_shell.number_measured_obs'))[0]
        if block.find_loop('_reflns_shell.percent_possible_all'):
            d['reflns_inner_percent_possible_obs'] = list(block.find_loop('_reflns_shell.percent_possible_all'))[0]
        if block.find_loop('_reflns_shell.pdbx_redundancy'):
            d['reflns_inner_pdbx_redundancy'] = list(block.find_loop('_reflns_shell.pdbx_redundancy'))[0]
        if block.find_loop('_reflns_shell.Rmerge_I_obs'):
            d['reflns_inner_pdbx_Rmerge_I_obs'] = list(block.find_loop('_reflns_shell.Rmerge_I_obs'))[0]
        if block.find_loop('_reflns_shell.meanI_over_sigI_obs'):
            d['reflns_inner_pdbx_netI_over_sigmaI'] = list(block.find_loop('_reflns_shell.meanI_over_sigI_obs'))[0]
        if block.find_loop('_reflns_shell.pdbx_Rrim_I_all'):
            d['reflns_inner_pdbx_pdbx_Rrim_I_all'] = list(block.find_loop('_reflns_shell.pdbx_Rrim_I_all'))[0]
        if block.find_loop('_reflns_shell.pdbx_CC_half'):
            d['reflns_inner_pdbx_CC_half'] = list(block.find_loop('_reflns_shell.pdbx_CC_half'))[0]
    return d
    

def get_highres_stats(block, d):
    if block.find_loop('_reflns_shell.pdbx_ordinal'):
        high = len(list(block.find_loop('_reflns_shell.pdbx_ordinal'))) - 1
        if block.find_loop('_reflns_shell.d_res_high'):
            d['reflns_outer_d_resolution_high'] = list(block.find_loop('_reflns_shell.d_res_high'))[high]
        if block.find_loop('_reflns_shell.d_res_low'):
            d['reflns_outer_d_resolution_low'] = list(block.find_loop('_reflns_shell.d_res_low'))[high]
        if block.find_loop('_reflns_shell.number_measured_obs'):
            d['reflns_outer_number_obs'] = list(block.find_loop('_reflns_shell.number_measured_obs'))[high]
        if block.find_loop('_reflns_shell.percent_possible_all'):
            d['reflns_outer_percent_possible_obs'] = list(block.find_loop('_reflns_shell.percent_possible_all'))[high]
        if block.find_loop('_reflns_shell.pdbx_redundancy'):
            d['reflns_outer_pdbx_redundancy'] = list(block.find_loop('_reflns_shell.pdbx_redundancy'))[high]
        if block.find_loop('_reflns_shell.Rmerge_I_obs'):
            d['reflns_outer_pdbx_Rmerge_I_obs'] = list(block.find_loop('_reflns_shell.Rmerge_I_obs'))[high]
        if block.find_loop('_reflns_shell.meanI_over_sigI_obs'):
            d['reflns_outer_pdbx_netI_over_sigmaI'] = list(block.find_loop('_reflns_shell.meanI_over_sigI_obs'))[high]
        if block.find_loop('_reflns_shell.pdbx_Rrim_I_all'):
            d['reflns_outer_pdbx_pdbx_Rrim_I_all'] = list(block.find_loop('_reflns_shell.pdbx_Rrim_I_all'))[high]
        if block.find_loop('_reflns_shell.pdbx_CC_half'):
            d['reflns_outer_pdbx_CC_half'] = list(block.find_loop('_reflns_shell.pdbx_CC_half'))[high]
    return d

def assign_dataset_outcome(logger, dal, d_proc):
    resolution = d_proc['reflns_d_resolution_high']
    logger.info('highest resolution for processed datasets is {0!s} A'.format(resolution))
    d = {}
    try:
        if resolution < 2.0:
            d['data_collection_outcome'] = "success - high resolution"
        elif resolution >= 2.0 and resolution <= 2.8:
            d['data_collection_outcome'] = "success - medium resolution"
        else:
            d['data_collection_outcome'] = "success - low resolution"
    except TypeError:
        d['data_collection_outcome'] = "unknown"
    logger.info('updating data_collection_outcome accordingly: {0!s}'.format(d['data_collection_outcome']))
    u = dal.xray_dataset_table.update().values(d).where(
        dal.xray_dataset_table.c.dataset_id == d_proc['dataset_id'])
    dal.connection.execute(u)

def update_processing_outcome(logger, dal, processing_id, processing_outcome):
    logger.info('updating processing_outcome')
    d = {}
    d['processing_outcome'] = processing_outcome
    u = dal.xray_processing_table.update().values(d).where(
        dal.xray_processing_table.c.processing_id == processing_id)
    dal.connection.execute(u)

def get_process_stats_from_mmcif_as_dict(logger, dal,ciffile, mtzfile, logfile, mounted_crystal_code, proposal, session, run, pipeline, projectDir):
    dataset_id = get_dataset_id(dal, mounted_crystal_code, proposal, session, run)

    d = {   'dataset_id':           dataset_id,
            'mounted_crystal_code': mounted_crystal_code,
            'automatic_processed':  True,
            'staraniso':            False,
            'processing_mtz_file':  mtzfile,
            'processing_cif_file':  ciffile,
            'processing_log_file':  logfile }

    mtz = gemmi.read_mtz_file(mtzfile)
    logger.info('ciffile: {0!s}'.format(os.path.realpath(ciffile)))

    try:
        doc = gemmi.cif.read_file(ciffile)
    except ValueError:
        logger.error('gemmi ValueError in processdb/get_process_stats_from_mmcif_as_dict reading {0!s}'.format(ciffile))
        return d

    d = get_cell_sym_info(mtz, d)

    for block in doc:
        d = get_software_info(logger, block, d)
        d = get_overall_stats(block, d)
        d = get_lowres_stats(block, d)
        d = get_highres_stats(block, d)
        if '_manual' in pipeline:
            break   # there are 4 blocks in Data_2_autoPROC_TRUNCATE_all.cif
                    # but script edits the automatically processed files so that only the header of the first
                    # block remains; note: xia2 has two blocks, the first being a summary, the second
                    # containing details for all resolution shells

#        break   # only interested in first block; xia2 has a second, somewhat redundant block
    if os.path.isfile(os.path.join(projectDir, '1-process', mounted_crystal_code, '{0!s}-{1!s}-{2!s}'.format(proposal, session, run), pipeline,
                       'unmerged_total.cif')):
        d = read_mrfana_cif(logger, d)
    logger.info('process_dict: {0!s}'.format(d))
    return d

def read_mrfana_cif(logger, d):
    if os.path.isfile('unmerged_total.cif'):
        logger.info('reading MRFANA CIF file')
        doc = gemmi.cif.read_file('unmerged_total.cif')
        for block in doc:
            found_20 = False
            found_15 = False
            found_10 = False
            if block.find_loop('_reflns_shell.pdbx_ordinal'):
                print(len(list(block.find_loop('_reflns_shell.d_res_high'))))
                print(len(list(block.find_loop('_reflns_shell.meanI_over_sigI_obs'))))
                for n, isig in enumerate(list(block.find_loop('_reflns_shell.meanI_over_sigI_obs'))):
                    if float(isig) < 2.0 and not found_20:
                        d['resolution_high_2_0_sigma'] = list(block.find_loop('_reflns_shell.d_res_high'))[n - 1]
                        found_20 = True
                    if float(isig) < 1.5 and not found_15:
                        d['resolution_high_1_5_sigma'] = list(block.find_loop('_reflns_shell.d_res_high'))[n - 1]
                        found_15 = True
                    if float(isig) < 1.0 and not found_10:
                        d['resolution_high_1_0_sigma'] = list(block.find_loop('_reflns_shell.d_res_high'))[n - 1]
                        found_10 = True
    return d


def insert_into_xray_processing_table(logger, dal, d):
    logger.info('saving xray_processing_table to database')
    try:
#        print(d)
        ins = dal.xray_processing_table.insert().values(d)
#        print(ins.compile(dialect=sqlite.dialect()))
        dal.connection.execute(ins)
    except sqlalchemy.exc.IntegrityError as e:
        if "UNIQUE constraint failed" in str(e):
            logger.warning('entry exists; skipping'.format(d['mounted_crystal_code']))
        else:
            logger.error(str(e))

def get_result_list_of_dicts(result):
    result_list = []
    for entry in result:
        value_dict = {}
        for n, key in enumerate(entry.keys()):
            if entry[n] == None:
                value_dict[key] = ''
            else:
                value_dict[key] = entry[n]
        result_list.append(value_dict)
    return result_list


def get_processing_results_for_sample(logger, dal, sample):
    logger.info('selecting all auto-processing results from database for {0!s}'.format(sample))
    q = select([dal.xray_processing_table.c.processing_id,
                dal.xray_processing_table.c.cell_volume,
                dal.xray_processing_table.c.sym_lattice_point_group,
                dal.xray_processing_table.c.sym_lattice,
                dal.xray_processing_table.c.sym_point_group,
                dal.xray_processing_table.c.reflns_d_resolution_high,
                dal.xray_processing_table.c.autoproc_pipeline,
                dal.xray_processing_table.c.automatic_processed,
                dal.xray_processing_table.c.staraniso,
                dal.xray_processing_table.c.data_reduction_software,
                dal.xray_processing_table.c.data_scaling_software,
                dal.xray_processing_table.c.reflns_inner_pdbx_Rmerge_I_obs,
                dal.xray_processing_table.c.processing_mtz_file,
                dal.xray_processing_table.c.processing_log_file,
                dal.xray_processing_table.c.processing_cif_file
                ]).where(dal.xray_processing_table.c.mounted_crystal_code == sample)
    rp = dal.connection.execute(q)
    r = rp.fetchall()
#    print(dir(dal.xray_processing_table.c.cell_volume))
#    print(dal.xray_processing_table.c.cell_volume.type)
#    sys.exit()
    result_list = get_result_list_of_dicts(r)
    for d in result_list:
        x = float(d['cell_volume'])
#        print(x)
#    sys.exit()
    return result_list


def unselected_autoprocessing_result(logger, dal, sample):
    logger.info('step 1: unselecting all auto-processing results for {0!s}'.format(sample))
    d = {}
    d['selected'] = False
    u = dal.xray_processing_table.update().values(d).where(
        dal.xray_processing_table.c.mounted_crystal_code == sample)
    dal.connection.execute(u)
#    unselect_datasets(logger, dal, sample)

def unselect_datasets(logger, dal, sample):
    logger.info('step 1a: unselecting all datasets for {0!s}'.format(sample))
    print(select(dal.xray_dataset_table))
    sys.exit()
    d = {}
    d['selected'] = False
    u = dal.xray_dataset_table.update().values(d).where(
        dal.xray_dataset_table.c.mounted_crystal_code == sample)
#    print("hallo")
#    print(u)
#    print(u.compile(dialect=sqlite.dialect()))
#    print('============================')
    dal.connection.execute(u)

def get_dataset_id_for_selected_autoprocessing_result(dal, mounted_crystal_code):
    q = select([dal.xray_processing_table.c.dataset_id]).where(and_(
                dal.xray_processing_table.c.mounted_crystal_code == mounted_crystal_code,
                dal.xray_processing_table.c.selected == True))
    rp = dal.connection.execute(q)
    r = rp.fetchall()
    idx = r[0][0]
    return idx

def set_selected_autoprocessing_result(logger, dal, sample, best):
    logger.info('step 2: set auto-processing results for {0!s}'.format(sample))
    d = {}
    d['selected'] = True
    # removed as constraint because is currently Null, but is "" in best dictionary
    # dal.xray_processing_table.c.data_scaling_software == best['data_scaling_software'],
    u = dal.xray_processing_table.update().values(d).where(and_(
        dal.xray_processing_table.c.mounted_crystal_code == sample,
        dal.xray_processing_table.c.data_reduction_software == best['data_reduction_software'],
        dal.xray_processing_table.c.autoproc_pipeline == best['autoproc_pipeline'],
        dal.xray_processing_table.c.automatic_processed == best['automatic_processed'],
        dal.xray_processing_table.c.staraniso == best['staraniso']))
    logger.info(d)
    logger.info(best)
    logger.info(u.compile(dialect=sqlite.dialect()))
    dal.connection.execute(u)
#    select_dataset(logger, dal, sample)

def select_dataset(logger, dal, sample):
    logger.info('step 2b: select dataset for auto-processing results of {0!s}'.format(sample))
    dataset_id = get_dataset_id_for_selected_autoprocessing_result(dal, sample)
    d = {}
    d['selected'] = True
    u = dal.xray_dataset_table.update().values(d).where(
        dal.xray_dataset_table.c.dataset_id == dataset_id)
    dal.connection.execute(u)

def get_dataset_id_of_last_dataset(dal, mounted_crystal_code):
    q = select([dal.xray_dataset_table.c.dataset_id,
                sqlalchemy.func.max(dal.xray_dataset_table.c.data_collection_date)]).where(and_(
                dal.xray_dataset_table.c.mounted_crystal_code == mounted_crystal_code,
                dal.xray_dataset_table.c.is_dataset == True))
    rp = dal.connection.execute(q)
    r = rp.fetchall()
    idx = r[0][0]
    return idx

def select_last_dataset(logger, dal, sample):
    logger.warning('no processing files, but select last collected dataset for summary query...')
    unselect_datasets(logger, dal, sample)
    dataset_id = get_dataset_id_of_last_dataset(dal, sample)
    d = {}
    d['selected'] = True
    u = dal.xray_dataset_table.update().values(d).where(
        dal.xray_dataset_table.c.dataset_id == dataset_id)
    dal.connection.execute(u)

def get_master_file_run_list(logger, dal, sample):
    logger.info('reading xray_dataset information for {0!s}'.format(sample))
#    q = select([dal.xray_dataset_table.c.h5_master_file,
#                dal.xray_dataset_table.c.run]).where(and_(
#                dal.xray_dataset_table.c.mounted_crystal_code == sample,
#                dal.xray_dataset_table.c.proposal == proposal,
#                dal.xray_dataset_table.c.session == session,
#                dal.xray_dataset_table.c.is_dataset == True))
#    q = select([dal.xray_dataset_table.c.h5_master_file,
#                dal.xray_dataset_table.c.run]).where(
#                dal.xray_dataset_table.c.mounted_crystal_code == sample)
    q = select([dal.xray_dataset_table.c.h5_master_file,
                dal.xray_dataset_table.c.run]).where(and_(
                dal.xray_dataset_table.c.mounted_crystal_code == sample,
                dal.xray_dataset_table.c.is_dataset == True))
    rp = dal.connection.execute(q)
    result = rp.fetchall()
    return result
