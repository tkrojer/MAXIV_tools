import getopt
import glob
import sys
import os
from datetime import datetime

from rdkit import Chem
from rdkit.Chem import AllChem

import pandas as pd
import json

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'lib')))
import processlib
import investigationlib


def add_project_data_to_cif(logger, block, data, name):
    investigationlib.add_pdbx_investigation(logger, block, data, name)
    investigationlib.add_fraghub_investigation_frag_screening(logger, block, data, name)
    investigationlib.add_pdbx_investigation_series(logger, block)
    investigationlib.add_pdbx_entity_poly(logger, block, data)
    investigationlib.add_pdbx_entity_poly_link_loop(logger, block, data)
    investigationlib.add_pdbx_contact_author_loop(logger, block, data)
    investigationlib.add_citation(logger, block)
    investigationlib.add_citation_author_loop(logger, block, data)
    investigationlib.add_pdbx_audit_support_loop(logger, block, data)


def find_all_entity_nonpoly_items(logger, projectDir, df, block):
    logger.info('step 1: finding all nonpoly items in model mmcif files')
    entity_nonpoly_dict = {}
    nonpoly_list = []
    for index, row in df.iterrows():
        sample_id, comp_id, smiles, pdb_id = investigationlib.get_row_items_from_df(logger, row)
        structure = investigationlib.get_model_mmcif(logger, projectDir, sample_id)
        if not structure:
            continue
        entity_nonpoly_dict = investigationlib.get_entity_nonpoly_from_structure_as_dict(
            logger, sample_id, structure, entity_nonpoly_dict)
    if entity_nonpoly_dict:
        nonpoly_list = investigationlib.prepare_nonpoly_list(logger, entity_nonpoly_dict)
        logger.info('found the following nonpoly items: {0!s}'.format(nonpoly_list))
        block = investigationlib.add_to_pdbx_entity_nonpoly_loop(logger, block, nonpoly_list)
    else:
        logger.warning('did not find any nonpoly items in model files; this is possible, but very unlikely!!!')
    logger.info('step 1: finished!')
    return block, nonpoly_list, entity_nonpoly_dict


def get_entity_nonpoly_link_list(logger, block, nonpoly_list, entity_nonpoly_dict):
    entity_nonpoly_link_list = []
    for sample_id in entity_nonpoly_dict:
        tmp = []
        for comp_id in entity_nonpoly_dict[sample_id]:
            for item in nonpoly_list:
                if comp_id == item['comp_id']:
                    tmp.append(item['entity_id'])
        tmp.sort()
        if tmp and tmp not in entity_nonpoly_link_list:
            entity_nonpoly_link_list.append(tmp)
    block = investigationlib.add_to_entity_nonpoly_link_loop(logger, block, entity_nonpoly_link_list)
    return entity_nonpoly_link_list


def find_entity_frag_library(logger, df, block):
    logger.info('step 2: finding fragment compound information')
    logger.info('')
    frag_library_list = []
    frag_library_entity_id = 1
    compound_list = []
    for index, row in df.iterrows():
        sample_id, comp_id, smiles, pdb_id = investigationlib.get_row_items_from_df(logger, row)
        print(smiles)
        if smiles:
            if comp_id not in compound_list:
                compound_list.append(comp_id)
                mol = investigationlib.get_rdkit_mol_from_smiles(smiles)
                frag_library_dict = investigationlib.get_frag_library_dict(
                    logger, mol, smiles, comp_id, frag_library_entity_id)
                frag_library_list.append(frag_library_dict)
                frag_library_entity_id += 1
        else:
            logger.warning('no CompoundID >{0!s}< and/or Smiles >{1!s}< associated with sample'.format(comp_id, smiles))
    if frag_library_list:
        logger.info('found {0!s} fragments'.format(len(frag_library_list)))
        block = investigationlib.add_to_fraghub_entity_frag_library_loop(logger, block, frag_library_list)
    else:
        logger.error('did not find information about any soaked fragments!')
    logger.info('step 2: finished')
    return block, frag_library_list


def get_fraghub_entity_frag_library_link(logger, block, frag_library_list):
    # this will be, at least for the time being the same as the entity_id in frag_library_list
    # because at the moment, we only record a single compound for each crystal
    # however, this may change in the future
    block = investigationlib.add_to_fraghub_entity_frag_library_link(logger, block, frag_library_list)
    return block


def get_pdbx_investigation_exp_nonpoly_dict(logger, entity_nonpoly_dict, entity_nonpoly_link_list, nonpoly_list):
    entity_nonpoly_link_dict = {}
    for sample_id in entity_nonpoly_dict:
        tmp = []
        for comp_id in entity_nonpoly_dict[sample_id]:
            for item in nonpoly_list:
                if comp_id == item['comp_id']:
                    tmp.append(item['entity_id'])
        tmp.sort()
        if tmp:
            entity_nonpoly_link_dict[sample_id] = str(entity_nonpoly_link_list.index(tmp)+1)
        else:
            entity_nonpoly_link_dict[sample_id] = '?'
    return entity_nonpoly_link_dict



def make_pdbx_investigation_exp_loop(logger, block, df, entity_nonpoly_link_dict, projectDir):
    logger.info('step 3: preparing _pdbx_investigation_exp loop')
    pdbx_investigation_exp_list = []
    for index, row in df.iterrows():
        sample_id, comp_id, smiles, pdb_id = investigationlib.get_row_items_from_df(logger, row)
        d = investigationlib.get_pdbx_investigation_exp_dict(index, sample_id, entity_nonpoly_link_dict)
        structure = investigationlib.get_model_mmcif(logger, projectDir, sample_id)
        if structure:
            d['exp_method'] = "X-RAY DIFFRACTION"
        if pdb_id:
            d['db'] = 'PDB'
            d['db_acc'] = pdb_id
        pdbx_investigation_exp_list.append(d)
    block = investigationlib.add_to_pdbx_investigation_exp_loop(logger, block, pdbx_investigation_exp_list)
    logger.info('step 3: finished')


def make_fraghub_investigation_frag_screening_exp_loop(logger, block, df, frag_library_list, projectDir):
    logger.info('step 4: preparing _fraghub_investigation_frag_screening_exp_loop loop')
    # here one can annotate each bound fragment, however, since we do not record this information
    # there will only be 3 categories (hit, miss - no fragment, miss - unknown (i.e. no model file exists)
    fraghub_investigation_frag_screening_exp_list = []
    counter = 1
    for index, row in df.iterrows():
        sample_id, comp_id, smiles, pdb_id = investigationlib.get_row_items_from_df(logger, row)
        structure = investigationlib.get_model_mmcif(logger, projectDir, sample_id)
        if pdb_id:
            instance = 1
            ligand_list = investigationlib.get_ligands_from_structure_as_list(structure)
            for ligand in ligand_list:
                d = investigationlib.get_fraghub_investigation_frag_screening_exp_dict(str(counter), index, frag_library_list, comp_id,
                                                                      str(instance))
                d['hit'] = "hit"
                d['hit_assessment'] = "refined"
                d['details'] = "On-site binding"
                fraghub_investigation_frag_screening_exp_list.append(d)
                instance += 1
                counter += 1

        elif structure:
            instance = 1
            d = investigationlib.get_fraghub_investigation_frag_screening_exp_dict(str(counter), index, frag_library_list, comp_id, str(instance))
            d['hit'] = "miss"
            d['hit_assessment'] = "manual"
            d['details'] = "Fragment unobserved"
            fraghub_investigation_frag_screening_exp_list.append(d)
            counter += 1
        else:
            instance = 1
            d = investigationlib.get_fraghub_investigation_frag_screening_exp_dict(str(counter), index, frag_library_list, comp_id, str(instance))
            d['hit'] = "miss"
            d['hit_assessment'] = "automatic"
            d['details'] = "unkown"
            fraghub_investigation_frag_screening_exp_list.append(d)
            counter += 1
    block = investigationlib.add_to_fraghub_investigation_frag_screening_exp_loop(logger, block, fraghub_investigation_frag_screening_exp_list)
    logger.info('step 4: finished')



def main(argv):
    projectDir = ''
    fragmaxcsv = ''
    projectjson = ''

    logger = processlib.init_logger('prepare_investigation_cif.log')
    processlib.start_logging(logger, 'prepare_investigation_cif.py')

    try:
        opts, args = getopt.getopt(argv,"p:f:j:h",["project_dir=", "fragmax_csv=", "json=", "help"])
    except getopt.GetoptError:
#        processlib.usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
#            processlib.usage()
            sys.exit(2)
        elif opt in ("-p", "--project_dir"):
            projectDir = os.path.abspath(arg)
        elif opt in ("-f", "--fragmax_csv"):
            fragmaxcsv = os.path.abspath(arg)
        elif opt in ("-j", "--json"):
            projectjson = os.path.abspath(arg)

    name = 'FRAGMAX1'
    json_file = open(projectjson)
    data = json.load(json_file)

    df = pd.read_csv(fragmaxcsv)
    df = df.where(pd.notnull(df), None)         # turns nan into None
    doc, block = investigationlib.init_investigation_cif(name)


    add_project_data_to_cif(logger, block, data, name)



    block, nonpoly_list, entity_nonpoly_dict = find_all_entity_nonpoly_items(logger, projectDir, df, block)
    entity_nonpoly_link_list = get_entity_nonpoly_link_list(logger, block, nonpoly_list, entity_nonpoly_dict)
    block, frag_library_list = find_entity_frag_library(logger, df, block)
    block = get_fraghub_entity_frag_library_link(logger, block, frag_library_list)
    entity_nonpoly_link_dict = get_pdbx_investigation_exp_nonpoly_dict(logger, entity_nonpoly_dict, entity_nonpoly_link_list, nonpoly_list)
    make_pdbx_investigation_exp_loop(logger, block, df, entity_nonpoly_link_dict, projectDir)
    make_fraghub_investigation_frag_screening_exp_loop(logger, block, df, frag_library_list, projectDir)
    doc.write_file('investigation.cif')
    logger.info('finished preparing investigation.cif file')

if __name__ == '__main__':
    main(sys.argv[1:])
