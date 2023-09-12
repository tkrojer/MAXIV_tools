from datetime import datetime
from sqlalchemy import (MetaData, Table, Column, Integer, Numeric, String,
        DateTime, ForeignKey, Boolean, create_engine, UniqueConstraint)

class DataAccessLayer:
    connection = None
    engine = None
    conn_string = None
    metadata = MetaData()

    """
    class generates all tables in database
    """

    project_table = Table('project_table',
        metadata,
        Column('project_id', Integer(), primary_key=True),
        Column('project_name', String(255), index=True),
        Column('proposal_number', String(50)),
        Column('project_description', String(255)),
        Column('project_directory', String(255)),
        Column('protein_acronym', String(50), unique=True),
        Column('protein_name', String(255)),
        Column('external_url', String(255)),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('protein_acronym', name='unique_protein_acronym')
    )

    author_table = Table('author_table',
        metadata,
        Column('author_id', Integer(), primary_key=True),
        Column('name_first', String(255)),
        Column('name_first_initials', String(255)),
        Column('name_last', String(255)),
        Column('address_1', String(255)),
        Column('address_2', String(255)),
        Column('address_3', String(255)),
        Column('email', String(255)),
        Column('role', String(255)),
        Column('identifier_ORCID', String(255)),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now)
    )

    funding_source_table = Table('funding_source_table',
        metadata,
        Column('funding_id', Integer(), primary_key=True),
        Column('funding_organization', String(255)),
        Column('country', String(255)),
        Column('grant_number', String(255)),
        Column('details', String(255)),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('funding_organization', 'grant_number', name='unique_grant')
    )

    protein_batch_table = Table('protein_batch_table',
        metadata,
        Column('protein_batch_id', Integer(), primary_key=True),
        Column('protein_batch_name', String(255)),
        Column('protein_acronym', String(50), ForeignKey('project_table.protein_acronym')),
        Column('protein_batch_supplier_name', String(255)),
        Column('protein_batch_uniprot_id', String(255)),
        Column('protein_batch_sequence', String(2550)),
        Column('protein_batch_expression_host', String(255)),
        Column('protein_batch_source_organism', String(255)),
        Column('protein_batch_vector', String(255)),
        Column('protein_batch_buffer', String(255)),
        Column('protein_batch_comp_id', String(255)),
        Column('protein_batch_concentration', Numeric(12, 2)),
        Column('protein_batch_concentration_unit_id', ForeignKey('unit_table.unit_id')),
        Column('protein_batch_date_received', DateTime(), default=datetime.now),
        Column('protein_batch_comment', String(2550)),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('protein_batch_name', name='unique_protein_batch_name')
    )

    crystal_screen_table = Table('crystal_screen_table',
        metadata,
        Column('crystal_screen_id', Integer(), primary_key=True),
        Column('crystal_screen_name', String(50), index=True),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('crystal_screen_name', name='unique_idx')
    )

    crystal_screen_condition_table = Table('crystal_screen_condition_table',
        metadata,
        Column('crystal_screen_condition_id', Integer(), primary_key=True),
        Column('crystal_screen_id', ForeignKey('crystal_screen_table.crystal_screen_id')),
        Column('crystal_screen_condition', String(255)),
        Column('crystal_screen_chem_comp_ids', String(255)),                # can be supplied comma separated
        Column('crystal_screen_chem_comp_concentrations', String(255)),     # in crystal_screen excel file
        Column('crystal_screen_chem_comp_units', String(255)),              #
        Column('crystal_screen_chem_comp_phs', String(255)),                #
        Column('crystal_screen_cocrystal_compound_code', String(255)),
        Column('crystal_screen_cocrystal_compound_concentration', String(255)),
        Column('crystal_screen_cocrystal_compound_concentration_unit_id', String(255)),
        Column('crystal_screen_well', String(12)),
        Column('crystal_screen_row', String(12)),
        Column('crystal_screen_column', String(12)),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('crystal_screen_id', 'crystal_screen_row', 'crystal_screen_column', name='unique_idx')
    )

    crystal_plate_table = Table('crystal_plate_table',
        metadata,
        Column('crystal_plate_id', Integer(), primary_key=True),
        Column('crystal_plate_barcode', String(255), index=True),
        Column('protein_batch_id', ForeignKey('protein_batch_table.protein_batch_id')),
        Column('protein_batch_concentration', Numeric(12, 2)),
        Column('protein_batch_concentration_unit_id', ForeignKey('unit_table.unit_id')),
        Column('protein_batch_buffer', String(255)),
        Column('protein_history', String(255)),
        Column('comment', String(255)),
        Column('crystal_screen_id', ForeignKey('crystal_screen_table.crystal_screen_id')),
        Column('method_id', ForeignKey('crystallization_method_table.method_id')),
        Column('temperature', Numeric(12, 2)),
        Column('plate_type_id', ForeignKey('plate_type_table.plate_type_id')),
        Column('reservoir_volume', Numeric(12, 2)),
        Column('reservoir_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('start_row', String(8)),
        Column('end_row', String(8)),
        Column('start_column', String(8)),
        Column('end_column', String(8)),
        Column('subwell_01_protein_volume', Numeric(12, 2)),
        Column('subwell_01_protein_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('subwell_01_reservoir_volume', Numeric(12, 2)),
        Column('subwell_01_reservoir_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('subwell_01_seed_volume', Numeric(12, 2)),
        Column('subwell_01_seed_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('subwell_00_protein_volume', Numeric(12, 2)),
        Column('subwell_00_protein_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('subwell_00_reservoir_volume', Numeric(12, 2)),
        Column('subwell_00_reservoir_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('subwell_00_seed_volume', Numeric(12, 2)),
        Column('subwell_00_seed_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('subwell_02_protein_volume', Numeric(12, 2)),
        Column('subwell_02_protein_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('subwell_02_reservoir_volume', Numeric(12, 2)),
        Column('subwell_02_reservoir_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('subwell_02_seed_volume', Numeric(12, 2)),
        Column('subwell_02_seed_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('subwell_03_protein_volume', Numeric(12, 2)),
        Column('subwell_03_protein_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('subwell_03_reservoir_volume', Numeric(12, 2)),
        Column('subwell_03_reservoir_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('subwell_03_seed_volume', Numeric(12, 2)),
        Column('subwell_03_seed_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('crystal_plate_barcode', name='unique_crystal_plate_barcode')
    )

    marked_crystals_table = Table('marked_crystals_table',
        metadata,
        Column('marked_crystal_id', Integer(), primary_key=True),
        Column('marked_crystal_code', String(255), index=True),
        Column('crystal_plate_id', ForeignKey('crystal_plate_table.crystal_plate_id')),
        Column('crystal_screen_condition_id', ForeignKey('crystal_screen_condition_table.crystal_screen_condition_id')),
        Column('crystal_plate_barcode', String(55)),
        Column('crystal_plate_row', String(55)),
        Column('crystal_plate_column', String(55)),
        Column('crystal_plate_subwell', String(55)),
        Column('crystal_plate_well', String(55)),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('marked_crystal_code', name='unique_marked_crystal_code')
    )

    soak_plate_table = Table('soak_plate_table',
        metadata,
        Column('soak_plate_id', Integer(), primary_key=True),
        Column('soak_plate_name', String(55)),
        Column('compound_batch_code', ForeignKey('compound_batch_table.compound_batch_code')),
        Column('compound_plate_name', String(55)),
        Column('soak_plate_type', String(20)),
        Column('soak_plate_type_id', ForeignKey('plate_type_table.plate_type_id')),
        Column('soak_plate_row', String(55)),
        Column('soak_plate_column', String(55)),
        Column('soak_plate_well', String(55)),
        Column('soak_plate_subwell', String(55)),
        Column('base_buffer', String(255)),
        Column('base_buffer_chem_comp_id', String(55)),
        Column('base_buffer_volume', Numeric(12, 2)),
        Column('base_buffer_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('compound_volume', Numeric(12, 2)),
        Column('compound_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('soak_plate_name', 'soak_plate_well', name='unique_name_well')
    )

    soaked_crystals_table = Table('soaked_crystals_table',
        metadata,
        Column('soaked_crystals_id', Integer(), primary_key=True),
        Column('soak_plate_name', String(55)),
        Column('soak_plate_id', ForeignKey('soak_plate_table.soak_plate_id')),
        Column('marked_crystal_code', String(255)),
        Column('marked_crystal_id', ForeignKey('marked_crystals_table.marked_crystal_id')),
        Column('compound_appearance', String(50)),
        Column('crystal_appearance', String(50)),
        Column('status', String(50)),   # soak_success, soak_fail, mount_success, mount_fail
        Column('soak_plate_type', String(20)),
        Column('soak_plate_well', String(55)),
        Column('soak_solution_volume', Numeric(12, 2)),
        Column('soak_solution_volume_unit', String(8)),
        Column('soak_solution_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('soak_method', String(55)),
        Column('soak_temperature', Numeric(12, 2)),
        Column('soak_temperature_unit_id', ForeignKey('unit_table.unit_id')),
        Column('comment', String(255)),
        Column('soak_datetime', DateTime(), default=datetime.now),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('soak_datetime', 'marked_crystal_id', name='unique_soak_datetime_crystal_id')
    )

    mounted_crystals_table = Table('mounted_crystals_table',
        metadata,
        Column('mounted_crystal_id', Integer(), primary_key=True),
        Column('mounted_crystal_code', String(255), index=True),
        Column('pin_barcode', String(55)),
        Column('puck_name', String(55)),
        Column('puck_position', Integer()),
        Column('mount_datetime', DateTime(), default=datetime.now),
        Column('mount_temperature', Numeric(12, 2)),
        Column('mount_temperature_unit_id', ForeignKey('unit_table.unit_id')),
        Column('shipment', String(55)),
        Column('marked_crystal_code', String(255)),
        Column('marked_crystal_id', ForeignKey('marked_crystals_table.marked_crystal_id')),
        Column('compound_appearance', String(50)),
        Column('crystal_appearance', String(50)),
        Column('cryo', String(55)),
        Column('cryo_chem_comp_code', String(55)),
        Column('cryo_volume', Numeric(12, 2)),
        Column('cryo_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('manual_mounted_crystal_code', String(55)),
        Column('data_collection_comment', String(55)),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('mount_datetime', name='unique_mount_datetime')
    )

    xray_dataset_table = Table('xray_dataset_table',
        metadata,
        Column('dataset_id', Integer(), primary_key=True),
        Column('mounted_crystal_id', ForeignKey('mounted_crystals_table.mounted_crystal_id')),
        Column('mounted_crystal_code', String(255)),
        Column('beamline', String(55)),
        Column('wavelength', String(55)),
        Column('proposal', String(55)),
        Column('session', String(55)),
        Column('run', String(55)),
        Column('h5_master_file', String(255)),
        Column('scicat_doi', String(255)),
        Column('dozor_plot', String(255)),
        Column('crystal_snapshot_1', String(255)),
        Column('crystal_snapshot_2', String(255)),
        Column('crystal_snapshot_3', String(255)),
        Column('crystal_snapshot_4', String(255)),
        Column('data_collection_outcome', String(255)),
        Column('data_collection_date', DateTime(), default=datetime.now),
        Column('data_collection_comment', String(255)),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('mounted_crystal_code', 'proposal', 'session', 'run', name='unique_xpsr')
    )

    xray_processing_table = Table('xray_processing_table',
        metadata,
        Column('processing_id', Integer(), primary_key=True),
        Column('dataset_id', ForeignKey('xray_dataset_table.dataset_id')),
        Column('mounted_crystal_code', String(255)),
        Column('data_reduction_software', String(55)),
        Column('data_reduction_software_version', String(55)),
        Column('data_scaling_software', String(55)),
        Column('data_scaling_software_version', String(55)),
        Column('autoproc_pipeline', String(55)),
        Column('autoproc_pipeline_version', String(55)),
        Column('automatic_processed', Boolean(), default=True),
        Column('cell_length_a', Numeric(12, 2)),
        Column('cell_length_b', Numeric(12, 2)),
        Column('cell_length_c', Numeric(12, 2)),
        Column('cell_angle_alpha', Numeric(12, 2)),
        Column('cell_angle_beta', Numeric(12, 2)),
        Column('cell_angle_gamma', Numeric(12, 2)),
        Column('cell_volume', Numeric(12, 2)),
        Column('sym_lattice', String(55)),
        Column('sym_point_group', String(55)),
        Column('sym_lattice_point_group', String(55)),
        Column('sym_space_group', String(55)),
        Column('sym_Int_Tables_number', Integer()),
        Column('reflns_d_resolution_high', Numeric(12, 2)),
        Column('reflns_d_resolution_low', Numeric(12, 2)),
        Column('resolution_high_1_0_sigma', Numeric(12, 2)),
        Column('resolution_high_1_5_sigma', Numeric(12, 2)),
        Column('resolution_high_2_0_sigma', Numeric(12, 2)),
        Column('reflns_number_obs', Integer()),
        Column('reflns_percent_possible_obs', Numeric(12, 2)),
        Column('reflns_pdbx_redundancy', Numeric(12, 2)),
        Column('reflns_pdbx_Rmerge_I_obs', Numeric(12, 2)),
        Column('reflns_pdbx_netI_over_sigmaI', Numeric(12, 2)),
        Column('reflns_pdbx_pdbx_Rrim_I_all', Numeric(12, 2)),
        Column('reflns_Rmeas_all', Numeric(12, 2)),
        Column('reflns_ISa_all', Numeric(12, 2)),
        Column('reflns_pdbx_CC_half', Numeric(12, 2)),
        Column('staraniso', Boolean(), default=False),
        Column('staraniso_version', String(55)),
        Column('reflns_B_iso_Wilson_estimate', Numeric(12, 2)),
        Column('reflns_pdbx_aniso_diffraction_limit_1', Numeric(12, 2)),
        Column('reflns_pdbx_aniso_diffraction_limit_2', Numeric(12, 2)),
        Column('reflns_pdbx_aniso_diffraction_limit_3', Numeric(12, 2)),
        Column('reflns_pdbx_percent_possible_ellipsoidal', Numeric(12, 2)),
        Column('reflns_pdbx_percent_possible_spherical', Numeric(12, 2)),
        Column('reflns_inner_d_resolution_high', Numeric(12, 2)),
        Column('reflns_inner_d_resolution_low', Numeric(12, 2)),
        Column('reflns_inner_number_obs', Integer()),
        Column('reflns_inner_percent_possible_obs', Numeric(12, 2)),
        Column('reflns_inner_pdbx_redundancy', Numeric(12, 2)),
        Column('reflns_inner_pdbx_Rmerge_I_obs', Numeric(12, 2)),
        Column('reflns_inner_pdbx_netI_over_sigmaI', Numeric(12, 2)),
        Column('reflns_inner_pdbx_pdbx_Rrim_I_all', Numeric(12, 2)),
        Column('reflns_inner_Rmeas_all', Numeric(12, 2)),
        Column('reflns_inner_pdbx_CC_half', Numeric(12, 2)),
        Column('reflns_outer_d_resolution_high', Numeric(12, 2)),
        Column('reflns_outer_d_resolution_low', Numeric(12, 2)),
        Column('reflns_outer_number_obs', Integer()),
        Column('reflns_outer_percent_possible_obs', Numeric(12, 2)),
        Column('reflns_outer_pdbx_redundancy', Numeric(12, 2)),
        Column('reflns_outer_pdbx_Rmerge_I_obs', Numeric(12, 2)),
        Column('reflns_outer_pdbx_netI_over_sigmaI', Numeric(12, 2)),
        Column('reflns_outer_pdbx_pdbx_Rrim_I_all', Numeric(12, 2)),
        Column('reflns_outer_Rmeas_all', Numeric(12, 2)),
        Column('reflns_outer_pdbx_CC_half', Numeric(12, 2)),
        Column('processing_directory_maxiv', String(255)),
        Column('processing_directory', String(255)),
        Column('processing_mtz_file', String(255)),
        Column('processing_log_file', String(255)),
        Column('processing_cif_file', String(255)),
        Column('processing_date', DateTime(), default=datetime.now),
        Column('processing_comment', String(255)),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('mounted_crystal_code',
                         'data_reduction_software',
                         'data_scaling_software',
                         'autoproc_pipeline',
                         'automatic_processed',
                         'staraniso',
                         name='unique_processing')
    )

    xray_initial_refinement_table = Table('xray_initial_refinement_table',
        metadata,
        Column('initial_refinement_id', Integer(), primary_key=True),
        Column('processing_id', ForeignKey('xray_processing_table.processing_id')),
        Column('mounted_crystal_code', String(255)),
        Column('refinement_software', String(55)),
        Column('initial_refinement_pipeline', String(55)),
        Column('input_pdb_file', String(255)),
        Column('input_cif_file', String(255)),
        Column('initial_refinement_directory', String(255)),
        Column('initial_refinement_mtz_file', String(255)),
        Column('initial_refinement_pdb_file', String(255)),
        Column('initial_refinement_cif_file', String(255)),
        Column('rfree_mtz_file', String(255)),
        Column('cell_length_a', Numeric(12, 2)),
        Column('cell_length_b', Numeric(12, 2)),
        Column('cell_length_c', Numeric(12, 2)),
        Column('cell_angle_alpha', Numeric(12, 2)),
        Column('cell_angle_beta', Numeric(12, 2)),
        Column('cell_angle_gamma', Numeric(12, 2)),
        Column('sym_space_group', String(55)),
        Column('sym_Int_Tables_number', Integer()),
        Column('refine_ls_d_res_high', Numeric(12, 2)),
        Column('refine_ls_d_res_low', Numeric(12, 2)),
        Column('refine_ls_R_factor_R_work', Numeric(12, 2)),
        Column('refine_ls_R_factor_R_free', Numeric(12, 2)),
        Column('refine_ls_percent_reflns_R_free', Numeric(12, 2)),
        Column('refine_r_bond_refined_d', Numeric(12, 2)),
        Column('refine_r_angle_refined_deg', Numeric(12, 2)),
        Column('initial_refinement_date', DateTime(), default=datetime.now),
        Column('initial_refinement_comment', String(255)),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('mounted_crystal_code',
                         'initial_refinement_pipeline',
                         name='unique_initial_refinement')
    )

    # xray_pandda_run_table
    # xray_pandda_analyse_event_table


    xray_refinement_table = Table('xray_refinement_table',
        metadata,
        Column('refinement_id', Integer(), primary_key=True),
        Column('initial_refinement_id', ForeignKey('xray_initial_refinement_table.initial_refinement_id')),
        Column('mounted_crystal_code', String(255)),
        Column('refinement_software', String(55)),
        Column('refinement_status', String(55)),
        Column('initial_refinement_date', DateTime(), default=datetime.now),
        Column('input_pdb_file', String(255)),
        Column('input_cif_file', String(255)),
        Column('refinement_directory', String(255)),
        Column('refinement_mtz_file', String(255)),
        Column('refinement_pdb_file', String(255)),
        Column('refinement_cif_file', String(255)),
        Column('rfree_mtz_file', String(255)),
        Column('cell_length_a', Numeric(12, 2)),
        Column('cell_length_b', Numeric(12, 2)),
        Column('cell_length_c', Numeric(12, 2)),
        Column('cell_angle_alpha', Numeric(12, 2)),
        Column('cell_angle_beta', Numeric(12, 2)),
        Column('cell_angle_gamma', Numeric(12, 2)),
        Column('sym_space_group', String(55)),
        Column('sym_Int_Tables_number', Integer()),
        Column('refine_ls_d_res_high', Numeric(12, 2)),
        Column('refine_ls_d_res_low', Numeric(12, 2)),
        Column('refine_ls_R_factor_R_work', Numeric(12, 2)),
        Column('refine_ls_R_factor_R_free', Numeric(12, 2)),
        Column('refine_ls_percent_reflns_R_free', Numeric(12, 2)),
        Column('refine_r_bond_refined_d', Numeric(12, 2)),
        Column('refine_r_angle_refined_deg', Numeric(12, 2)),
        Column('refinement_date', DateTime(), default=datetime.now),
        Column('refinement_comment', String(255)),
        Column('refinement_cycle', Integer()),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('mounted_crystal_code', name='unique_refinement')
    )


    xray_model_table = Table('xray_model_table',
        metadata,
        Column('model_id', Integer(), primary_key=True),
        Column('refinement_id', ForeignKey('xray_refinement_table.refinement_id')),
        Column('mounted_crystal_code', String(255)),
        Column('model_status', String(55)),
        Column('model_directory', String(255)),
        Column('model_mtz_file', String(255)),
        Column('model_pdb_file', String(255)),
        Column('model_cif_file', String(255)),
        Column('model_comment', String(255)),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('mounted_crystal_code', name='unique_refinement')
    )

    # xray_model_observations_table

    pdb_table = Table('pdb_table',
        metadata,
        Column('pdb_table_id', Integer(), primary_key=True),
        Column('model_id', ForeignKey('xray_model_table.model_id')),
        Column('mounted_crystal_code', String(255)),
        Column('pdb_code', String(55)),
        Column('group_code', String(55)),
        Column('archived_data_url', String(255)),
        Column('external_data_url', String(255)),
        Column('external_doi', String(255)),
        Column('group_deposition', Boolean(), default=False),
        Column('deposition_date', DateTime(), default=datetime.now),
        Column('deposition_comment', String(255)),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('pdb_code', name='unique_pdb_code')
    )

    compound_table = Table('compound_table',
        metadata,
        Column('compound_id', Integer(), primary_key=True),
        Column('compound_code', String(50), index=True, unique=True),
        Column('chemical_name', String(255)),
        Column('smiles', String(255)),
        Column('inchi', String(255)),
        Column('cas', String(255)),
        Column('vendor', String(55)),
        Column('vendor_id', String(55)),
        Column('chem_comp_code', String(55)),
        Column('formula', String(255)),
        Column('formula_weight', Numeric(12, 2)),
        Column('plate_name', String(55)),
        Column('cocktail', Boolean(), default=False),
        Column('covalent', Boolean(), default=False),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('compound_code', name='unique_compound_code')
    )

    compound_batch_table = Table('compound_batch_table',
        metadata,
        Column('compound_batch_id', Integer(), primary_key=True),
        Column('compound_batch_code', String(50), index=True),
        Column('compound_code', ForeignKey('compound_table.compound_code')),
        Column('library_name', String(55)),
        Column('compound_plate_name', String(55)),
        Column('compound_plate_row', String(8)),
        Column('compound_plate_column', String(8)),
        Column('compound_plate_well', String(8)),
        Column('compound_plate_type', ForeignKey('plate_type_table.plate_type_id')),
        Column('solvent', ForeignKey('compound_table.compound_code')),
        Column('solvent_chem_comp_code', String(55)),
        Column('compound_concentration', Numeric(12, 2)),
        Column('compound_concentration_unit', String(50)),
        Column('compound_concentration_unit_id', ForeignKey('unit_table.unit_id')),
        Column('compound_volume', Numeric(12, 2)),
        Column('compound_volume_unit', String(50)),
        Column('compound_volume_unit_id', ForeignKey('unit_table.unit_id')),
        Column('comment', String(255)),
        Column('compound_batch_active', Boolean(), default=False),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('compound_batch_code', name='unique_compound_batch_code')
    )

    wwpdb_chem_comp_table = Table('wwpdb_chem_comp_table',
        metadata,
        Column('chem_comp_id', Integer(), primary_key=True),
        Column('chem_comp_code', String(50), index=True),
        Column('name', String(255)),
        Column('popular_name', String(255)),
        Column('type', String(255)),
        Column('pdbx_type', String(255)),
        Column('formula', String(255)),
        Column('formula_weight', Numeric(12, 2)),
        Column('smiles', String(255)),
        Column('inchi', String(255)),
        Column('created_on', DateTime(), default=datetime.now),
        Column('updated_on', DateTime(), default=datetime.now, onupdate=datetime.now),
        UniqueConstraint('chem_comp_code', name='unique_chem_comp_code')
    )

    unit_table = Table('unit_table',
        metadata,
        Column('unit_id', Integer(), primary_key=True),
        Column('unit_symbol', String(50), index=True),
        Column('unit_name', String(50)),
        Column('conversion_factor', Numeric(12, 12)),
        Column('type', String(50)),
        UniqueConstraint('unit_symbol', name='unique_unit_symbol')
    )

    plate_type_table = Table('plate_type_table',
        metadata,
        Column('plate_type_id', Integer(), primary_key=True),
        Column('plate_type_name', String(20), index=True),
        UniqueConstraint('plate_type_name', name='unique_plate_type_name')
    )

    crystallization_method_table = Table('crystallization_method_table',
        metadata,
        Column('method_id', Integer(), primary_key=True),
        Column('method_name', String(50), index=True),
        UniqueConstraint('method_name', name='unique_method_name')
    )

    soak_method_table = Table('soak_method_table',
        metadata,
        Column('method_id', Integer(), primary_key=True),
        Column('method_name', String(50), index=True),
        UniqueConstraint('method_name', name='unique_method_name')
    )

    space_group_table = Table('space_group_table',
        metadata,
        Column('space_group_id', Integer(), primary_key=True),
        Column('space_group_name', String(50), index=True),
        Column('space_group_number', Integer()),
        UniqueConstraint('space_group_name', name='unique_space_group_name')
    )

    gene_src_table = Table('gene_src_table',
        metadata,
        Column('pdbx_gene_src_id', Integer(), primary_key=True),
        Column('pdbx_gene_src_scientific_name', String(50), index=True),
        Column('pdbx_gene_src_ncbi_taxonomy_id', Integer()),
        UniqueConstraint('pdbx_gene_src_scientific_name', name='unique_pdbx_gene_src_scientific_name')
    )

    wwpdb_allowed_values_table = Table('wwpdb_allowed_values_table',
        metadata,
        Column('wwpdb_allowed_values_id', Integer(), primary_key=True),
        Column('pdbx_category', String(50), index=True),
        Column('allowed_value', String(255)),
        UniqueConstraint('pdbx_category', 'allowed_value', name='unique_pdbx_gene_src_scientific_name')
    )

    def db_init(self, conn_string):
        self.engine = create_engine('sqlite:///' + conn_string or self.conn_string)
        self.metadata.create_all(self.engine)
        self.connection = self.engine.connect()

dal = DataAccessLayer()
