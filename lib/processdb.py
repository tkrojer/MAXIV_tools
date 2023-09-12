import gemmi

import sqlalchemy
from sqlalchemy.sql import select
from sqlalchemy import and_

def get_mounted_crystal_id(dal, sample):
    q = select([dal.mounted_crystals_table.c.mounted_crystal_id]).where(
        dal.mounted_crystals_table.c.mounted_crystal_code == sample)
    rp = dal.connection.execute(q)
    result = rp.fetchall()
    return result[0][0]

def get_d_xray_dataset_table_dict(logger, dal, sample, proposal, session, beamline, run, create_date, master, dozor_plot):
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

    return d_xray_dataset_table_dict

def insert_into_xray_dataset_table(logger, dal, d):
    logger.info('trying to insert into xray_dataset_table')
    try:
        ins = dal.xray_dataset_table.insert().values(d)
        dal.connection.execute(ins)
    except sqlalchemy.exc.IntegrityError as e:
        if "UNIQUE constraint failed" in str(e):
            logger.warning('entry exists (time soaked {0!s}); skipping'.format(d['mounted_crystal_code']))
        else:
            logger.error(str(e))

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

def get_highres_stats(block, d):
    if block.find_loop('_reflns_shell.pdbx_ordinal'):
        high = len(list(block.find_loop('_reflns_shell.pdbx_ordinal'))) - 1
        if block.find_loop('_reflns_shell.d_res_high'):
            d['reflns_inner_d_resolution_high'] = list(block.find_loop('_reflns_shell.d_res_high'))[high]
        if block.find_loop('_reflns_shell.d_res_low'):
            d['reflns_inner_d_resolution_low'] = list(block.find_loop('_reflns_shell.d_res_low'))[high]
        if block.find_loop('_reflns_shell.number_measured_obs'):
            d['reflns_inner_number_obs'] = list(block.find_loop('_reflns_shell.number_measured_obs'))[high]
        if block.find_loop('_reflns_shell.percent_possible_all'):
            d['reflns_inner_percent_possible_obs'] = list(block.find_loop('_reflns_shell.percent_possible_all'))[high]
        if block.find_loop('_reflns_shell.pdbx_redundancy'):
            d['reflns_inner_pdbx_redundancy'] = list(block.find_loop('_reflns_shell.pdbx_redundancy'))[high]
        if block.find_loop('_reflns_shell.Rmerge_I_obs'):
            d['reflns_inner_pdbx_Rmerge_I_obs'] = list(block.find_loop('_reflns_shell.Rmerge_I_obs'))[high]
        if block.find_loop('_reflns_shell.meanI_over_sigI_obs'):
            d['reflns_inner_pdbx_netI_over_sigmaI'] = list(block.find_loop('_reflns_shell.meanI_over_sigI_obs'))[high]
        if block.find_loop('_reflns_shell.pdbx_Rrim_I_all'):
            d['reflns_inner_pdbx_pdbx_Rrim_I_all'] = list(block.find_loop('_reflns_shell.pdbx_Rrim_I_all'))[high]
        if block.find_loop('_reflns_shell.pdbx_CC_half'):
            d['reflns_inner_pdbx_CC_half'] = list(block.find_loop('_reflns_shell.pdbx_CC_half'))[high]
    return d
    

def get_lowres_stats(block, d):
    if block.find_loop('_reflns_shell.pdbx_ordinal'):
        if block.find_loop('_reflns_shell.d_res_high'):
            d['reflns_outer_d_resolution_high'] = list(block.find_loop('_reflns_shell.d_res_high'))[0]
        if block.find_loop('_reflns_shell.d_res_low'):
            d['reflns_outer_d_resolution_low'] = list(block.find_loop('_reflns_shell.d_res_low'))[0]
        if block.find_loop('_reflns_shell.number_measured_obs'):
            d['reflns_outer_number_obs'] = list(block.find_loop('_reflns_shell.number_measured_obs'))[0]
        if block.find_loop('_reflns_shell.percent_possible_all'):
            d['reflns_outer_percent_possible_obs'] = list(block.find_loop('_reflns_shell.percent_possible_all'))[0]
        if block.find_loop('_reflns_shell.pdbx_redundancy'):
            d['reflns_outer_pdbx_redundancy'] = list(block.find_loop('_reflns_shell.pdbx_redundancy'))[0]
        if block.find_loop('_reflns_shell.Rmerge_I_obs'):
            d['reflns_outer_pdbx_Rmerge_I_obs'] = list(block.find_loop('_reflns_shell.Rmerge_I_obs'))[0]
        if block.find_loop('_reflns_shell.meanI_over_sigI_obs'):
            d['reflns_outer_pdbx_netI_over_sigmaI'] = list(block.find_loop('_reflns_shell.meanI_over_sigI_obs'))[0]
        if block.find_loop('_reflns_shell.pdbx_Rrim_I_all'):
            d['reflns_outer_pdbx_pdbx_Rrim_I_all'] = list(block.find_loop('_reflns_shell.pdbx_Rrim_I_all'))[0]
        if block.find_loop('_reflns_shell.pdbx_CC_half'):
            d['reflns_outer_pdbx_CC_half'] = list(block.find_loop('_reflns_shell.pdbx_CC_half'))[0]
    return d

def get_process_stats_from_mmcif_as_dict(logger, dal,ciffile, mtzfile, logfile, mounted_crystal_code, proposal, session, run):

    dataset_id = get_dataset_id(dal, mounted_crystal_code, proposal, session, run)

    d = {   'dataset_id':           dataset_id,
            'mounted_crystal_code': mounted_crystal_code,
            'automatic_processed':  True,
            'staraniso':            False,
            'processing_mtz_file':  mtzfile,
            'processing_cif_file':  ciffile,
            'processing_log_file':  logfile }

    mtz = gemmi.read_mtz_file(mtzfile)
    doc = gemmi.cif.read_file(ciffile)

    d = get_cell_sym_info(mtz, d)

    for block in doc:
        d = get_software_info(logger, block, d)
        d = get_overall_stats(block, d)
        d = get_lowres_stats(block, d)
        d = get_highres_stats(block, d)

    return d


def insert_into_xray_processing_table(logger, dal, d):
    logger.info('saving xray_processing_table to database')
    try:
        ins = dal.xray_processing_table.insert().values(d)
        dal.connection.execute(ins)
    except sqlalchemy.exc.IntegrityError as e:
        if "UNIQUE constraint failed" in str(e):
            logger.warning('entry exists (time soaked {0!s}); skipping'.format(d['soak_datetime']))
        else:
            logger.error(str(e))
