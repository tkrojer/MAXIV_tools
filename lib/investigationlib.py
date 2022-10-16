import os
import sys

from rdkit import Chem
from rdkit.Chem import Descriptors

import gemmi

import requests


def init_investigation_cif(name):
    doc = gemmi.cif.Document()
    doc.add_new_block(name)
    block = doc.sole_block()
    return doc, block


def add_pdbx_investigation(logger, block, data, name):
    logger.info('adding _pdbx_investigation items')
    title = '?'
    details = '?'
    for item in data:
        if item == 'DepositionTitle':
            title = data[item]
        if item == 'Description':
            details = data[item]
    block.set_pair('_pdbx_investigation.id', name)
    block.set_pair('_pdbx_investigation.title', gemmi.cif.quote(title))
    block.set_pair('_pdbx_investigation.details', gemmi.cif.quote(details))
    block.set_pair('_pdbx_investigation.external_url', gemmi.cif.quote('?'))


def add_fraghub_investigation_frag_screening(logger, block, data, name):
    logger.info('adding _fraghub_investigation_frag_screening items')
    uniprot_id = '?'
    for item in data:
        if item == 'Entities':
            for _list in data[item]:
                uniprot_id = _list['UniprotID']
                break
    block.set_pair('_fraghub_investigation_frag_screening.investigation_id', name)
    block.set_pair('_fraghub_investigation_frag_screening.campaign_id', '1')
    block.set_pair('_fraghub_investigation_frag_screening.target', 'MID2')
    block.set_pair('_fraghub_investigation_frag_screening.common_unp', uniprot_id)
    block.set_pair('_fraghub_investigation_frag_screening.facility', gemmi.cif.quote('MAX IV'))
    block.set_pair('_fraghub_investigation_frag_screening.proc_pipeline', 'FRAGMAXAPP')
    block.set_pair('_fraghub_investigation_frag_screening.internal_id', '20211214')


def add_pdbx_investigation_series(logger, block):
    logger.info('adding _pdbx_investigation_series items')
    block.set_pair('_pdbx_investigation_series.frag_screening_campaign_id', '1')
    block.set_pair('_pdbx_investigation_series.id', '1')
    block.set_pair('_pdbx_investigation_series.fragment_lib', 'FRAGMAXLIB')
    block.set_pair('_pdbx_investigation_series.fragment_batch', 'v1.0')
    block.set_pair('_pdbx_investigation_series.details', '?')


def init_pdbx_entity_poly_loop(block):
    loop = block.init_loop('_pdbx_entity_poly.', [
                'entity_id',
                'type',
                'src_method',
                'seq_one_letter_code',
                'db',
                'db_acc'
        ])
    return block, loop


def add_pdbx_entity_poly(logger, block, data):
    logger.info('adding _pdbx_entity_poly loop items')
    block, loop = init_pdbx_entity_poly_loop(block)
    for item in data:
        if item == 'Entities':
            for n, _list in enumerate(data[item]):
                uniprot_id = _list['UniprotID']
                sequence = _list['Sequence']
                loop.add_row([
                    str(n+1),
                    gemmi.cif.quote("polypeptide(L)"),
                    "man",
                    sequence,
                    "UNP",
                    uniprot_id
                    ])


def init_pdbx_entity_poly_link_loop(block):
    loop = block.init_loop('_pdbx_entity_poly_link.', [
                'id',
                'entity_id'
        ])
    return block, loop


def add_pdbx_entity_poly_link_loop(logger, block, data):
    logger.info('adding _pdbx_entity_poly_link loop items')
    block, loop = init_pdbx_entity_poly_link_loop(block)
    for item in data:
        if item == 'Entities':
            for n, _list in enumerate(data[item]):
                loop.add_row([
                    '1',
                    str(n+1)
                ])


def init_pdbx_contact_author_loop(block):
    loop = block.init_loop('_pdbx_contact_author.', [
                'id',
                'name_first',
                'name_last',
                'address_1',
                'address_2',
                'address_3',
                'email',
                'role'
        ])
    return block, loop


def add_pdbx_contact_author_loop(logger, block, data):
    logger.info('adding _pdbx_contact_author loop items')
    block, loop = init_pdbx_contact_author_loop(block)
    for item in data:
        if item == 'Scientists':
            for n, _list in enumerate(data[item]):
                loop.add_row([
                    str(n+1),
                    gemmi.cif.quote(_list['FirstName']),
                    gemmi.cif.quote(_list['LastName']),
                    gemmi.cif.quote(_list['OrganizationName']),
                    gemmi.cif.quote(_list['Street']),
                    gemmi.cif.quote(_list['ZIPCode'] + ', ' + _list['City'] + ', ' +_list['Country']),
                    gemmi.cif.quote(_list['Email']),
                    gemmi.cif.quote(_list['Role'])
                    ])


def add_citation(logger, block):
    logger.info('adding _citation items')
    block.set_pair('_citation.abstract', '?')
    block.set_pair('_citation.abstract_id_CAS', '?')
    block.set_pair('_citation.book_id_ISBN', '?')
    block.set_pair('_citation.book_publisher', '?')
    block.set_pair('_citation.book_publisher_city', '?')
    block.set_pair('_citation.book_title', '?')
    block.set_pair('_citation.coordinate_linkage', '?')
    block.set_pair('_citation.country', '?')
    block.set_pair('_citation.database_id_Medline', '?')
    block.set_pair('_citation.details', '?')
    block.set_pair('_citation.id', 'primary')
    block.set_pair('_citation.journal_abbrev', '?')
    block.set_pair('_citation.journal_id_ASTM', '?')
    block.set_pair('_citation.journal_id_CSD', '?')
    block.set_pair('_citation.journal_id_ISSN', '?')
    block.set_pair('_citation.journal_full', '?')
    block.set_pair('_citation.journal_issue', '?')
    block.set_pair('_citation.journal_volume', '?')
    block.set_pair('_citation.language', '?')
    block.set_pair('_citation.page_first', '?')
    block.set_pair('_citation.page_last', '?')
    block.set_pair('_citation.title', '?')
    block.set_pair('_citation.year', '?')
    block.set_pair('_citation.pdbx_database_id_DOI', '?')
    block.set_pair('_citation.pdbx_database_id_PubMed', '?')


def init_citation_author_loop(block):
    loop = block.init_loop('_citation_author.', [
                'citation_id',
                'name',
                'ordinal',
                'identifier_ORCID'
        ])
    return block, loop


def add_citation_author_loop(logger, block, data):
    logger.info('adding _citation_author loop items')
    block, loop = init_citation_author_loop(block)
    for item in data:
        if item == 'Authors':
            for n, _list in enumerate(data[item]):
                loop.add_row([
                    'primary',
                    gemmi.cif.quote(_list['Name']),
                    str(n+1),
                    gemmi.cif.quote(_list['ORCID'])
                ])


def init_pdbx_audit_support_loop(block):
    loop = block.init_loop('_citation_author.', [
                'ordinal',
                'funding_organization',
                'country',
                'grant_number',
                'details'
        ])
    return block, loop


def add_pdbx_audit_support_loop(logger, block, data):
    logger.info('adding _pdbx_audit_support loop items')
    block, loop = init_pdbx_audit_support_loop(block)
    for item in data:
        if item == 'Funding':
            for n, _list in enumerate(data[item]):
                loop.add_row([
                    str(n+1),
                    gemmi.cif.quote(_list['Organization']),
                    gemmi.cif.quote(_list['Country']),
                    gemmi.cif.quote(_list['GrantNumber']),
                    '?'
                    ])


def get_row_items_from_df(logger, row):
    sample_id = row['SampleID']
    comp_id = row['CompoundID']
    smiles = row['Smiles']
    pdb_id = row['PDBID']
    logger.info('current sample ' + sample_id)
    return sample_id, comp_id, smiles, pdb_id


def get_model_mmcif(logger, projectDir, sample_id):
    mmcif = None
#    print(os.path.join(projectDir, sample_id, sample_id + '_model.cif'))
    if os.path.isfile(os.path.join(projectDir, sample_id, sample_id + '_model.cif')):
        logger.info('found model mmcif file')
        mmcif = os.path.join(projectDir, sample_id, sample_id + '_model.cif')
    else:
        logger.error('cannot find model mmcif for sample')
    if mmcif:
        structure = gemmi.read_structure(mmcif)
    else:
        structure = None
    return structure


def get_entity_nonpoly_from_structure_as_dict(logger, sample_id, structure, entity_nonpoly_dict):
    ligand_residue_names = ['LIG', 'DRG', 'XXX', 'YYY', '220', '241', '188', '234', '250']
    entity_nonpoly_dict[sample_id] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if not gemmi.find_tabulated_residue(residue.name).is_amino_acid():
                    if residue.name not in ligand_residue_names and \
                        residue.name not in entity_nonpoly_dict[sample_id]:
                        logger.info('found entity_nonpoly item ' + residue.name)
                        entity_nonpoly_dict[sample_id].append(residue.name)
    return entity_nonpoly_dict


def get_ligands_from_structure_as_list(structure):
    ligand_residue_names = ['LIG', 'DRG', 'XXX', 'YYY', '220', '241', '188', '234', '250']
    ligand_list = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.name in ligand_residue_names:
                    ligand_list.append(residue.name)
    return ligand_list


def get_empty_nonpoly_dict(n, comp_id):
    d = {
        'entity_id': str(n+1),
        'name': '?',
        'comp_id': comp_id,
        'formula': '?',
        'formula_weight': '?',
        'inchi_descriptor': '?',
        'cas_identifier': '?',
    }
    return d


def prepare_nonpoly_list(logger, entity_nonpoly_dict):
    comp_list = []
    for sample in entity_nonpoly_dict:
        for comp_id in entity_nonpoly_dict[sample]:
            if comp_id not in comp_list:
                comp_list.append(comp_id)
    nonpoly_list = []
    for n, comp_id in enumerate(comp_list):
        d = get_empty_nonpoly_dict(n, comp_id)
        letter = str(comp_id)[0].lower()
        cif_file = '/Applications/ccp4-7.1//lib/data/monomers/{0!s}/{1!s}.cif'.format(letter, str(comp_id).upper())
        if os.path.isfile(cif_file):
            doc = gemmi.cif.read_file(cif_file)
            for block in doc:
                if block.find_loop('_chem_comp.name'):
                    d['name'] = str(list(block.find_loop('_chem_comp.name'))[0]).replace("'", '').rstrip()
                if block.find_loop('_pdbx_chem_comp_descriptor.descriptor'):
                    for item in list(block.find_loop('_pdbx_chem_comp_descriptor.descriptor')):
                        if item.replace('"', '').startswith('InChI'):
                            d['inchi_descriptor'] = item.replace('"', '')
                            mol = Chem.MolFromInchi(d['inchi_descriptor'])
                            d['formula'] = formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                            d['formula_weight'] = str(round(Descriptors.MolWt(mol), 1))
        nonpoly_list.append(d)
    return nonpoly_list


def init_pdbx_entity_nonpoly_loop(block):
    loop = block.init_loop('_pdbx_entity_nonpoly.', [   'entity_id',
                                                        'name',
                                                        'comp_id',
                                                        'formula',
                                                        'formula_weight',
                                                        'inchi_descriptor',
                                                        'cas_identifier'
        ])
    return block, loop

def add_to_pdbx_entity_nonpoly_loop(logger, block, nonpoly_list):
    logger.info('adding _pdbx_entity_nonpoly to investigation cif file...')
    block, loop = init_pdbx_entity_nonpoly_loop(block)
    for d in nonpoly_list:
        logger.info('-> ' + d['name'])
        loop.add_row([  d['entity_id'],
                        gemmi.cif.quote(d['name']),
                        d['comp_id'],
                        gemmi.cif.quote(d['formula']),
                        d['formula_weight'],
                        gemmi.cif.quote(d['inchi_descriptor']),
                        d['cas_identifier']
                  ])
    return block


def get_rdkit_mol_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol


def get_empty_frag_library_dict(comp_id, frag_library_entity_id):
    ligand_residue_names = ['LIG', 'DRG', 'XXX', 'YYY', '220', '241', '188', '234', '250']
    if comp_id == 'VT00220':
        lig_id = '220'
    elif comp_id == 'VT00241':
        lig_id = '241'
    elif comp_id == 'VT00259':
        lig_id = 'XXX'
    elif comp_id == 'VT00188':
        lig_id = '188'
    elif comp_id == 'VT00234':
        lig_id = '234'
    elif comp_id == 'VT00250':
        lig_id = '250'
    elif comp_id == 'VT00258':
        lig_id = 'XXX'
    elif comp_id == 'VT00259':
        lig_id = 'XXX'
    elif comp_id == 'VT00431':
        lig_id = 'XXX'
    else:
        lig_id = 'LIG'


    d = {
        'entity_id':        str(frag_library_entity_id),
        'parent_id':        '0',
        'series_id':        '1',
        'name':             comp_id,
        'details':          '?',
        'comp_id':          lig_id,
        'formula':          '',
        'formula_weight':   '',
        'inchi_descriptor': '',
        'cas_identifier':   '?'
    }
    return d


def get_frag_library_dict(logger, mol, smiles, comp_id, frag_library_entity_id):
    d = get_empty_frag_library_dict(comp_id, frag_library_entity_id)
    d = get_properties(logger, mol, smiles, comp_id, d)
    return d


def get_properties(logger, mol, smiles, comp_id, d):
    d['formula'] = Chem.rdMolDescriptors.CalcMolFormula(mol)
    d['formula_weight'] = str(round(Descriptors.MolWt(mol), 1))
    d['inchi_descriptor'] = Chem.MolToInchi(mol)
#    iupac = get_iupac_name_from_smiles(smiles)
    logger.info('converting {0!s} smiles {1!s} to {2!s}'.format(comp_id, smiles, d['formula_weight']))
    logger.info('formula: {0!s} - MW: {1!s}'.format(d['formula'], d['formula_weight']))
    return d


def get_iupac_name_from_smiles(smiles):
    CACTUS = "https://cactus.nci.nih.gov/chemical/structure/{0}/{1}"
    rep = "iupac_name"
    url = CACTUS.format(smiles, rep)
    try:
        response = requests.get(url)
        response.raise_for_status()
        iupac = response.text
    except requests.exceptions.HTTPError:
        iupac = '?'
    return iupac


def init_fraghub_entity_frag_library_loop(block):
    loop = block.init_loop('_fraghub_entity_frag_library.', [   'entity_id',
                                                                'parent_id',
                                                                'series_id',
                                                                'name',
                                                                'details',
                                                                'comp_id',
                                                                'formula',
                                                                'formula_weight',
                                                                'inchi_descriptor',
                                                                'cas_identifier'  ])
    return block, loop

def add_to_fraghub_entity_frag_library_loop(logger, block, frag_library_list):
    logger.info('adding _fraghub_entity_frag_library to investigation cif file...')
    block, loop = init_fraghub_entity_frag_library_loop(block)
    for d in frag_library_list:
        logger.info('-> ' + d['name'])
        loop.add_row([  d['entity_id'],
                        d['parent_id'],
                        d['series_id'],
                        gemmi.cif.quote(d['name']),
                        d['details'],
                        d['comp_id'],
                        gemmi.cif.quote(d['formula']),
                        d['formula_weight'],
                        gemmi.cif.quote(d['inchi_descriptor']),
                        d['cas_identifier']
                  ])
    return block


def init_entity_nonpoly_link_loop(block):
    loop = block.init_loop('_pdbx_entity_nonpoly_link.', [  'id',
                                                            'entity_id'  ])
    return block, loop

def add_to_entity_nonpoly_link_loop(logger, block, entity_nonpoly_link_list):
    logger.info('adding _entity_nonpoly_link to investigation cif file...')
    block, loop = init_entity_nonpoly_link_loop(block)
    for n, item in enumerate(entity_nonpoly_link_list):
        for i in item:
            loop.add_row([  str(n+1),
                            i
                  ])
    return block


def init_fraghub_entity_frag_library_link(block):
    loop = block.init_loop('_fraghub_entity_frag_library_link.', [  'id',
                                                                    'entity_id'  ])
    return block, loop

def add_to_fraghub_entity_frag_library_link(logger, block, frag_library_list):
    # see comment in prepare_investigation_cif.py
    # at the moment there will be only a single compound for each crystal
    logger.info('adding _fraghub_entity_frag_library_link to investigation cif file...')
    block, loop = init_fraghub_entity_frag_library_link(block)
    for item in frag_library_list:
        loop.add_row([  item['entity_id'],
                        item['entity_id']
                  ])
    return block


def get_pdbx_investigation_exp_dict(index, sample_id, entity_nonpoly_link_dict):
    if sample_id in entity_nonpoly_link_dict:
        nonpoly_link = entity_nonpoly_link_dict[sample_id]
    else:
        nonpoly_link = '?'

    d = {
            'id': str(index + 1),
            'series_id': '1',
            'exp_acc': sample_id,
            'exp_poly': '1',
            'exp_nonpoly': nonpoly_link,
            'exp_method': '?',
            'db': '?',
            'db_acc': '?',
            'exp_details': '?',
            'external_url': '?'
    }
    return d


def init_pdbx_investigation_exp_loop(block):
    loop = block.init_loop('_pdbx_investigation_exp.', [
                'id',
                'series_id',
                'exp_acc',
                'exp_poly',
                'exp_nonpoly',
                'exp_method',
                'db',
                'db_acc',
                'exp_details',
                'external_url'
        ])
    return block, loop


def add_to_pdbx_investigation_exp_loop(logger, block, pdbx_investigation_exp_list):
    logger.info('adding _pdbx_investigation_exp to investigation cif file...')
    block, loop = init_pdbx_investigation_exp_loop(block)
    for d in pdbx_investigation_exp_list:
        loop.add_row([
            d['id'],
            d['series_id'],
            d['exp_acc'],
            d['exp_poly'],
            d['exp_nonpoly'],
            gemmi.cif.quote(d['exp_method']),
            d['db'],
            d['db_acc'],
            d['exp_details'],
            d['external_url']
        ])
    return block


def get_fraghub_investigation_frag_screening_exp_dict(counter, index, frag_library_list, comp_id, instance):
    # note:
    # exp_fragment and bound_fragments are identical since the script only deals with one fragment per crystal
    exp_fragment = '?'
    bound_fragments = '?'

    for d in frag_library_list:
        if d['name'] == comp_id:
            exp_fragment = d['entity_id']
            bound_fragments = d['entity_id']
            break

    if comp_id in frag_library_list:
        exp_fragment = frag_library_list
    d = {
        'id':               counter,
        'exp_id':           str(index + 1),
        'instance_id':      instance,
        'exp_fragment':     exp_fragment,
        'bound_fragments':  bound_fragments,
        'hit':              '',
        'hit_assessment':   '',
        'details':          ''
    }
    return d


def init_fraghub_investigation_frag_screening_exp_loop(block):
    loop = block.init_loop('_fraghub_investigation_frag_screening_exp.', [
        'id',
        'exp_id',
        'instance_id',
        'exp_fragment',
        'bound_fragments',
        'hit',
        'hit_assessment',
        'details'
    ])
    return block, loop


def add_to_fraghub_investigation_frag_screening_exp_loop(logger, block, fraghub_investigation_frag_screening_exp_list):
    logger.info('adding _fraghub_investigation_frag_screening_exp to investigation cif file...')
    block, loop = init_fraghub_investigation_frag_screening_exp_loop(block)
    for d in fraghub_investigation_frag_screening_exp_list:
        loop.add_row([
            d['id'],
            d['exp_id'],
            d['instance_id'],
            d['exp_fragment'],
            d['bound_fragments'],
            gemmi.cif.quote(d['hit']),
            gemmi.cif.quote(d['hit_assessment']),
            gemmi.cif.quote(d['details'])
        ])
    return block
