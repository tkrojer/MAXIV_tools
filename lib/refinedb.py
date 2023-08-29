import gemmi

import sqlalchemy
from sqlalchemy.sql import select
from sqlalchemy import and_


def get_processing_id(dal, data):
    mounted_crystal_code = data['mounted_crystal_code']
    run = data['run']
    proposal = data['proposal']
    session = data['session']
    autoproc_pipeline = data['pipeline']

    j = dal.xray_dataset_table.join(
            dal.xray_processing_table, dal.xray_dataset_table.c.dataset_id ==
                                       dal.xray_processing_table.c.dataset_id, isouter=True)

    q = select([dal.xray_processing_table.c.processing_id]).where(and_(
                dal.xray_dataset_table.c.proposal == proposal,
                dal.xray_dataset_table.c.session == session,
                dal.xray_dataset_table.c.run == run,
                dal.xray_processing_table.c.mounted_crystal_code == mounted_crystal_code,
                dal.xray_processing_table.c.autoproc_pipeline == autoproc_pipeline))

    q = q.select_from(j)
    rp = dal.connection.execute(q)
    r = rp.fetchall()
    idx = r[0][0]
    return idx

def get_xray_initial_refinement_table_dict(logger, dal, sample, data):
    logger.info('getting d_xray_initial_refinement_table_dict for {0!s}'.format(sample))
    process_id = get_processing_id(dal, data)

    d = {   'processing_id':        process_id,
            'mounted_crystal_code': sample }





    return d


def insert_into_xray_initial_refinement_table(logger, dal, d):
    logger.info('saving xxray_initial_refinement_table to database')
    try:
        ins = dal.xray_initial_refinement_table.insert().values(d)
        dal.connection.execute(ins)
    except sqlalchemy.exc.IntegrityError as e:
        if "UNIQUE constraint failed" in str(e):
            logger.warning('entry exists (time soaked {0!s}); skipping'.format(d['soak_datetime']))
        else:
            logger.error(str(e))
