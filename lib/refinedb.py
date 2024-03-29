import gemmi

import sqlalchemy
from sqlalchemy.sql import select
from sqlalchemy import and_

def get_processing_id(logger, dal, sample):
    q = select([dal.xray_processing_table.c.processing_id]).where(and_(
                dal.xray_processing_table.c.mounted_crystal_code == sample,
                dal.xray_processing_table.c.selected == True))
    rp = dal.connection.execute(q)
    result = rp.fetchall()
    logger.info('processing_id for {0!s} is {1!s}'.format(sample, result[0][0]))
    return result[0][0]

def get_d_xray_initial_refinement_table_dict(logger, dal, sample, software, pdbref, mtzin, cifref, mtzref):
    logger.info('getting xray_initial_refinement_table_dict for {0!s}'.format(sample))
    processing_id = get_processing_id(logger, dal, sample)
    if software == 'dimple':
        refinement_software = 'refmac'
    elif software == 'pipedream':
        refinement_software = 'buster'
    else:
        refinement_software = 'unknown'

    d = {
        'processing_id': processing_id,
        'mounted_crystal_code': sample,
        'initial_refinement_pipeline': software,
        'refinement_software': refinement_software,
        'input_pdb_file': pdbref,
        'input_mtz_file': mtzin,
        'input_cif_file': cifref,
        'input_mtz_free_file': mtzref
    }
    return d

def insert_update_xray_initial_refinement_table(logger, dal, d, sample, software):
    logger.info('saving xray_initial_refinement_table to database')
    try:
        ins = dal.xray_initial_refinement_table.insert().values(d)
        dal.connection.execute(ins)
    except sqlalchemy.exc.IntegrityError as e:
        if "UNIQUE constraint failed" in str(e):
            d.pop('mounted_crystal_code', 'mounted_crystal_code not found in dictionary')
            d.pop('initial_refinement_pipeline', 'initial_refinement_pipeline not found in dictionary')
            logger.warning('entry exists: {0!s}/; updating...'.format(sample, software))
            u = dal.xray_initial_refinement_table.update().values(d).where(and_(
                dal.xray_initial_refinement_table.c.mounted_crystal_code == sample,
                dal.xray_initial_refinement_table.c.initial_refinement_pipeline == software))
            dal.connection.execute(u)
        else:
            logger.error(str(e))

def unselected_initial_refinement_pipeline(logger, dal, sample):
    logger.info('step 1: unselecting all initial refinement results for {0!s}'.format(sample))
    d = {}
    d['selected'] = False
    u = dal.xray_initial_refinement_table.update().values(d).where(
        dal.xray_initial_refinement_table.c.mounted_crystal_code == sample)
    dal.connection.execute(u)


def set_selected_initial_refinement_pipeline(logger, dal, sample, software):
    logger.info('step 2: set initial refinement results for {0!s}'.format(sample))
    d = {}
    d['selected'] = True
    u = dal.xray_initial_refinement_table.update().values(d).where(and_(
        dal.xray_initial_refinement_table.c.mounted_crystal_code == sample,
        dal.xray_initial_refinement_table.c.initial_refinement_pipeline == software))
    dal.connection.execute(u)

def update_initial_refinement_outcome(logger, dal, sample, software, initial_refinement_outcome):
    logger.warning('setting outcome for {0!s} - {1!s} to {2!s}'.format(sample, software, initial_refinement_outcome))
    d = {}
    d['initial_refinement_outcome'] = initial_refinement_outcome
    u = dal.xray_initial_refinement_table.update().values(d).where(and_(
        dal.xray_initial_refinement_table.c.mounted_crystal_code == sample,
        dal.xray_initial_refinement_table.c.initial_refinement_pipeline == software))
    dal.connection.execute(u)



