import sys
import os
import getopt
import xlsxwriter
import logging
from datetime import datetime
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'lib')))
import ciftools

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

def get_data_from_cif(logger, model_cif, process_cif, validate_xml):
    d = ciftools.get_empty_dict()
    doc = ciftools.get_cif_as_doc(logger, model_cif)
    d = ciftools.get_refine_stats_from_doc_as_dict(logger, doc, d)
    d = ciftools.get_atom_count_baverage_as_dict(logger, model_cif, d)
    if process_cif:
        doc = ciftools.get_cif_as_doc(logger, process_cif)
    d = ciftools.get_process_stats_from_doc_as_dict(logger, doc, d)
    if validate_xml:
        d = ciftools.get_info_from_validation_xml_as_dict(logger, validate_xml, d)
    return d

def make_workbook(logger, model_cif, process_cif, validate_xml):
    d = get_data_from_cif(logger, model_cif, process_cif, validate_xml)
    # Create a new Excel file and add a worksheet
    simplified_workbook = xlsxwriter.Workbook('table_one.xlsx')
    simplified_worksheet = simplified_workbook.add_worksheet()
    print(d)
    # Data as observed from the simplified spreadsheet, adjusted to remove the middle column
    structure_data = [
        ("Data Collection", None),
        ("PDB ID", "XXX"),
        ("Beamline", "BioMAX"),
        ("Wavelength (Å)", "XXX"),
        ("Space Group", f"{d['sym_space_group']}"),
        ("Cell dimensions", None),
        ("a, b, c (Å)", f"{d['cell_length_a']}, {d['cell_length_b']}, {d['cell_length_c']}"),
        ("⍺, β, ɣ (degrees)", f"{d['cell_angle_alpha']}, {d['cell_angle_beta']}, {d['cell_angle_gamma']}"),
        ("Resolution range (Å)", f"{d['reflns_d_resolution_low']} - {d['reflns_d_resolution_high']} ({d['reflns_outer_d_resolution_low']} - {d['reflns_d_resolution_high']})"),
        ("No. measured intensities", f"{d['reflns_number_obs']}"),
        ("No. unique reflections", f"{d['reflns_number_obs']}"),
        ("Multiplicity", f"{d['reflns_pdbx_redundancy']} ({d['reflns_outer_pdbx_redundancy']})"),
        ("Mean I/σ(I)", f"{d['reflns_pdbx_netI_over_sigmaI']} ({d['reflns_outer_pdbx_netI_over_sigmaI']})"),
        ("Completeness (%)", f"{d['reflns_percent_possible_obs']} ({d['reflns_outer_percent_possible_obs']})"),
        ("Rmerge", f"{d['reflns_pdbx_Rmerge_I_obs']} ({d['reflns_outer_pdbx_Rmerge_I_obs']})"),
        ("CC1⁄2", f"{d['reflns_pdbx_CC_half']} ({d['reflns_outer_pdbx_CC_half']})"),
        ("Wilson B value (Å2)", f"{d['wilson_b_estimate']}"),
        ("Refinement", None),
        ("Resolution range (Å)", f"{d['refine_ls_d_res_low']} - {d['refine_ls_d_res_high']}"),
        ("Reflections: working/free", None),
        ("Rwork (%)", f"{d['refine_ls_R_factor_R_work']}"),
        ("Rfree (%)", f"{d['refine_ls_R_factor_R_free']}"),
        ("No. Atoms", None),
        ("Protein", f"{d['n_protein_atom']}"),
        ("Water", f"{d['n_water_atom']}"),
        ("Ligand", f"{d['n_ligand_atom']}"),
        ("Other", f"{d['n_other_atom']}"),
        ("B-factors", None),
        ("Protein", f"{d['bfac_protein']}"),
        ("Water", f"{d['bfac_water']}"),
        ("Ligand", f"{d['bfac_ligand']}"),
        ("Other", f"{d['bfac_other']}"),
        ("R. m. s. deviations", None),
        ("Bond lengths (Å)", f"{d['refine_r_bond_refined_d']}"),
        ("Bond angles (°)", f"{d['refine_r_angle_refined_deg']}"),
        ("Ramachandran plot (%)", None),
        ("favoured", f"{d['percent_rama_favoured']}"),
        ("outliers", f"{d['percent_rama_outliers']}")
    ]

    # Populate the worksheet with the structure and data without the middle column
    for row_num, (structure, value) in enumerate(structure_data):
        simplified_worksheet.write(row_num, 0, structure)  # Writing the structure/category
        if value is not None:
            simplified_worksheet.write(row_num, 1, value)  # Writing the value directly next to the label

    # Closing the workbook to save the file
    simplified_workbook.close()

def main(argv):
    now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    logger = init_logger(f"logfile_{now}.log")

    model_cif = None
    process_cif = None
    validate_xml = None

    try:
        opts, args = getopt.getopt(argv,"m:c:x:h",["model=", "cif=", "xml=", "help"])
    except getopt.GetoptError:
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            logger.info("ccp4-python table_one.py -m <model.mmcif> -c <process.cif> -x <validate.xml>")
            sys.exit(2)
        elif opt in ("-m", "--model"):
            model_cif = os.path.abspath(arg)
        elif opt in ("-c", "--cif"):
            process_cif = os.path.abspath(arg)
        elif opt in ("-x", "--xml"):
            validate_xml = os.path.abspath(arg)

    if model_cif:
        make_workbook(logger, model_cif, process_cif, validate_xml)


if __name__ == "__main__":
    main(sys.argv[1:])