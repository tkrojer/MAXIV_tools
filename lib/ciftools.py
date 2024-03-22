import gemmi
from xml.etree import ElementTree

def get_cif_as_doc(logger, ciffile):
    doc = None
    try:
        logger.info(f"reading {ciffile} as doc")
        doc = gemmi.cif.read_file(ciffile)
    except ValueError:
        logger.error(f'gemmi ValueError when reading {ciffile}')
    return doc

def get_process_stats_from_doc_as_dict(logger, doc, d):
    logger.info("getting data processing statistics from cif doc object")
    for block in doc:
#        d = get_software_info(logger, block, d)
        d = get_overall_stats(logger, block, d)
        d = get_lowres_stats(logger, block, d)
        d = get_highres_stats(logger, block, d)
    return d

def get_overall_stats(logger, block, d):
    if block.find_pair('_reflns.d_resolution_low'):
        logger.info("found overall information")
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

def get_lowres_stats(logger, block, d):
    if block.find_loop('_reflns_shell.pdbx_ordinal'):
        logger.info("found information for low resolution shell")
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

def get_highres_stats(logger, block, d):
    if block.find_loop('_reflns_shell.pdbx_ordinal'):
        logger.info("found information for high resolution shell")
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

def get_residue_categories():
    aa_list = ["GLY", "ALA", "VAL", "LEU", "ILE", "THR", "SER", "MET",
               "MSE", "CYS", "PRO", "PHE", "TYR", "TRP", "HIS", "LYS",
               "ARG", "ASP", "GLU", "ASN", "GLN"]

    hoh_list = ['HOH']

    lig_list = ['LIG', 'DRG']
    return aa_list, hoh_list, lig_list

def get_atom_count_baverage_as_dict(logger, model_mmcif, d):
    logger.info(f"reading {model_mmcif}")
#    structure = gemmi.read_structure(model_mmcif, merge_chain_parts=True)
    try:
        structure = gemmi.read_structure(model_mmcif)
    except ValueError:
        cif_block = gemmi.cif.read(model_mmcif)[0]
        structure = gemmi.make_structure_from_block(cif_block)

    aa_list, hoh_list, lig_list = get_residue_categories()

    d['n_protein_atom'] = 0
    d['n_water_atom'] = 0
    d['n_ligand_atom'] = 0
    d['n_other_atom'] = 0

    d['bfac_protein'] = 0
    d['bfac_water'] = 0
    d['bfac_ligand'] = 0
    d['bfac_other'] = 0

    bfac_list_protein = []
    bfac_list_water = []
    bfac_list_ligand = []
    bfac_list_other = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if residue.name in aa_list:
                        bfac_list_protein.append(atom.b_iso)
                        d['n_protein_atom'] += 1
                    elif residue.name in lig_list:
                        bfac_list_ligand.append(atom.b_iso)
                        d['n_ligand_atom'] += 1
                    elif residue.name in hoh_list:
                        bfac_list_water.append(atom.b_iso)
                        d['n_water_atom'] += 1
                    else:
                        bfac_list_other.append(atom.b_iso)
                        d['n_other_atom'] += 1

    if bfac_list_protein:
        d['bfac_protein'] = round((sum(bfac_list_protein) / len(bfac_list_protein)), 2)
    if bfac_list_water:
        d['bfac_water'] = round((sum(bfac_list_water) / len(bfac_list_water)), 2)
    if bfac_list_ligand:
        d['bfac_ligand'] = round((sum(bfac_list_ligand) / len(bfac_list_ligand)), 2)
    if bfac_list_other:
        d['bfac_other'] = round((sum(bfac_list_other) / len(bfac_list_other)), 2)

    return d

def get_refine_stats_from_doc_as_dict(logger, doc, d):
    logger.info("getting refinement statistics from cif doc object")
    for block in doc:
        if block.find_pair('_software.name'):
            d['refinement_program'] = block.find_pair('_software.name')[1]
        if block.find_pair('_refine.ls_R_factor_R_work'):
            d['refine_ls_R_factor_R_work'] = round(float(block.find_pair('_refine.ls_R_factor_R_work')[1]), 2)
        if block.find_pair('_refine.ls_R_factor_R_free'):
            d['refine_ls_R_factor_R_free'] = round(float(block.find_pair('_refine.ls_R_factor_R_free')[1]), 2)
        if block.find_pair('_refine.ls_d_res_high'):
            d['refine_ls_d_res_high'] = round(float(block.find_pair('_refine.ls_d_res_high')[1]), 2)
        if block.find_pair('_refine.ls_d_res_low'):
            d['refine_ls_d_res_low'] = round(float(block.find_pair('_refine.ls_d_res_low')[1]), 2)
        if block.find_pair('_cell.length_a'):
            d['cell_length_a'] = round(float(block.find_pair('_cell.length_a')[1]), 2)
        if block.find_pair('_cell.length_b'):
            d['cell_length_b'] = round(float(block.find_pair('_cell.length_b')[1]), 2)
        if block.find_pair('_cell.length_c'):
            d['cell_length_c'] = round(float(block.find_pair('_cell.length_c')[1]), 2)
        if block.find_pair('_cell.angle_alpha'):
            d['cell_angle_alpha'] = round(float(block.find_pair('_cell.angle_alpha')[1]), 2)
        if block.find_pair('_cell.angle_beta'):
            d['cell_angle_beta'] = round(float(block.find_pair('_cell.angle_beta')[1]), 2)
        if block.find_pair('_cell.angle_gamma'):
            d['cell_angle_gamma'] = round(float(block.find_pair('_cell.angle_gamma')[1]), 2)
        if block.find_pair('_symmetry.space_group_name_H-M'):
            d['sym_space_group'] = block.find_pair('_symmetry.space_group_name_H-M')[1]
        if block.find_pair('_symmetry.Int_Tables_number'):
            d['sym_Int_Tables_number'] = block.find_pair('_symmetry.Int_Tables_number')[1]
        if block.find_loop('_refine_ls_restr.type') and 'refinement_program' in d:
            if d['refinement_program'] == 'refmac':
                table = block.find('_refine_ls_restr.', ['type', 'dev_ideal'])
                d['refine_r_bond_refined_d'] = round(float(list(table.find_row('r_bond_refined_d'))[1]), 3)
                d['refine_r_angle_refined_deg'] = round(float(list(table.find_row('r_angle_refined_deg'))[1]), 3)
            if d['refinement_program'] == 'buster':
                table = block.find('_refine_ls_restr.', ['type', 'dev_ideal'])
                d['refine_r_bond_refined_d'] = round(float(list(table.find_row('t_bond_d'))[1]), 3)
                d['refine_r_angle_refined_deg'] = round(float(list(table.find_row('t_angle_deg'))[1]), 3)
    return d

def get_info_from_validation_xml_as_dict(logger, xml, d):
    logger.info(f"reading {xml} from onedep validation")
    tree = ElementTree.parse(xml)
    for item in tree.getroot():
        if item.tag == 'Entry':
            d['percent_rama_outliers'] = dict(item.items())['percent_rama_outliers']
            d['percent_rama_favoured'] = 100 - float(d['percent_rama_outliers'])
            d['wilson_b_estimate'] = round(float(dict(item.items())['WilsonBestimate']), 1)
    return d

def get_empty_dict():
    d = {
#        "Wavelength (Å)": "",
        'sym_space_group': "",
        'cell_length_a': "",
        'cell_length_b': "",
        'cell_length_c': "",
        'cell_angle_alpha': "",
        'cell_angle_beta': "",
        'cell_angle_gamma': "",
        'reflns_d_resolution_low': "",
        'reflns_d_resolution_high': "",
        'reflns_outer_d_resolution_low': "",
        'reflns_number_obs': "",
        "No. unique reflections": "",
        'reflns_pdbx_redundancy': "",
        'reflns_outer_pdbx_redundancy': "",
        'reflns_pdbx_netI_over_sigmaI': "",
        'reflns_outer_pdbx_netI_over_sigmaI': "",
        'reflns_percent_possible_obs': "",
        'reflns_outer_percent_possible_obs': "",
        'reflns_pdbx_Rmerge_I_obs': "",
        'reflns_outer_pdbx_Rmerge_I_obs': "",
        'reflns_pdbx_CC_half': "",
        'reflns_outer_pdbx_CC_half': "",
        'wilson_b_estimate': "",
        'refine_ls_d_res_low': "",
        'refine_ls_d_res_high': "",
        'refine_ls_R_factor_R_work': "",
        'refine_ls_R_factor_R_free': "",
        'n_protein_atom': "",
        'n_water_atom': "",
        'n_ligand_atom': "",
        'n_other_atom': "",
        'bfac_protein': "",
        'bfac_water': "",
        'bfac_ligand': "",
        'bfac_other': "",
        'refine_r_bond_refined_d': "",
        'refine_r_angle_refined_deg': "",
        'percent_rama_favoured': "",
        'percent_rama_outliers': ""
    }

    return d