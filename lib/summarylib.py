import xlsxwriter
from datetime import datetime
import json
import os
import sys
import glob
from PIL import Image

import seaborn as sns
import matplotlib.pyplot as plt
sns.set_theme(style="ticks", color_codes=True)
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'lib')))
import processlib
import refinelib


def init_workbook():
    now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    workbook = xlsxwriter.Workbook('summary_{0!s}.xlsx'.format(now))
    return workbook


def get_n_samples(fragmaxcsv):
    n_samples = 0
    for n,line in enumerate(open(fragmaxcsv)):
        sample = line.split(',')[0]
    n_samples = n + 1
    return n_samples


def n_autoproc_results(projectDir):
    n_autoproc_results = 0
    for cif in glob.glob(os.path.join(projectDir, '1-process', '*', '*', '*', 'process.cif')):
        n_autoproc_results += 1
    return n_autoproc_results


def get_blank_cif():
    cif = {
        'collection_date': '',
        'proposal': '',
        'session': '',
        'reso_low': '',
        'reso_high': '',
        'status': '',
        'pipeline': '',
        'percent_possible_obs': '',
        'Rmerge_I_obs_low': '',
        'meanI_over_sigI_obs_high': '',
        'unitcell': '',
        'space_group': '',
        'initref_software': '',
        'initref_spacegroup': '',
        'ls_d_res_high': '',
        'ls_R_factor_R_work': '',
        'ls_R_factor_R_free': '',
        'r_bond_refined_d': '',
        'r_angle_refined_deg': '',
        'blobs': ''
    }
    return cif

def add_cif_columns_to_df(df):
    df['collection_date'] = ' '
    df['proposal'] = ' '
    df['session'] = ' '
    df['reso_low'] = ' '
    df['reso_high'] = ' '
    df['status'] = ' '
    df['pipeline'] = ' '
    df['percent_possible_obs'] = ' '
    df['Rmerge_I_obs_low'] = ' '
    df['meanI_over_sigI_obs_high'] = ' '
#    df['unitcell'] = ' '
    df['space_group'] = ' '
    df['initref_software'] = ' '
    df['initref_spacegroup'] = ' '
    df['ls_d_res_high'] = ' '
    df['ls_R_factor_R_work'] = ' '
    df['ls_R_factor_R_free'] = ' '
    df['r_bond_refined_d'] = ' '
    df['r_angle_refined_deg'] = ' '
    df['blobs'] = ' '
    df['smiles'] = ' '
    df['sample_appearance'] = ' '
    df['cpd_appearance'] = ' '
    return df

def get_semi_blank_cif():
    cif = {
        'reso_low': '',
        'reso_high': '',
        'pipeline': '',
        'percent_possible_obs': '',
        'Rmerge_I_obs_low': '',
        'meanI_over_sigI_obs_high': '',
        'unitcell': '',
        'space_group': '',
        'initref_software': '',
        'initref_spacegroup': '',
        'ls_d_res_high': '',
        'ls_R_factor_R_work': '',
        'ls_R_factor_R_free': '',
        'r_bond_refined_d': '',
        'r_angle_refined_deg': '',
        'blobs': ''
    }
    return cif


def get_blank_init_refine_cif():
    cif = {
        'initref_software': '',
        'initref_spacegroup': '',
        'ls_d_res_high': '',
        'ls_R_factor_R_work': '',
        'ls_R_factor_R_free': '',
        'r_bond_refined_d': '',
        'r_angle_refined_deg': '',
        'blobs': ''
    }
    return cif


def get_json_as_dict(jsonfile):
    with open(jsonfile) as d:
        dictData = json.load(d)
    return dictData


def number_to_column(n):
    columnDict = {
        '1': 'A',
        '2': 'B',
        '3': 'C',
        '4': 'D',
        '5': 'E',
        '6': 'F',
        '7': 'G',
        '8': 'H',
        '9': 'I',
        '10': 'J',
        '11': 'K',
        '12': 'L',
        '13': 'M',
        '14': 'N',
        '15': 'O',
        '16': 'P',
        '17': 'Q',
        '18': 'R',
        '19': 'S',
        '20': 'T',
        '21': 'U',
        '22': 'V',
        '23': 'W',
        '24': 'X',
        '25': 'Y',
        '26': 'Z'
    }
    return columnDict[str(n)]


def get_point_group_dict():
    pgDict = {}
    return pgDict


def update_point_group_dict(cif, pgDict):
    pg = str(cif['lattice']) + str(cif['point_group'])
    if pg not in pgDict:
        pgDict[pg] = 0
    pgDict[pg] += 1
    return pgDict


def get_point_group_info_dict():
    pginfoDict = {
        'pointgroup': [],
        'unitcell_volume': [],
        'Rmerge_I_obs_low': []
    }
    return pginfoDict


def update_point_group_info_dict(cif, pginfoDict):
    pg = str(cif['lattice']) + str(cif['point_group'])
    ucv = int(cif['unitcell_volume'])
    rmergelow = float(cif['Rmerge_I_obs_low'])
    pginfoDict['pointgroup'].append(pg)
    pginfoDict['unitcell_volume'].append(ucv)
    pginfoDict['Rmerge_I_obs_low'].append(rmergelow)
    return pginfoDict


#def get_point_group_rmerge_low_dict():
#    pgrmergelowDict = {
#        'pointgroup': [],
#        'Rmerge_I_obs_low': []
#    }
#    return pgrmergelowDict


#def update_point_group_rmerge_low_dict(cif, pgrmergelowDict):
#    pg = cif['lattice'] + cif['point_group']
#    rmergelow = float(cif['Rmerge_I_obs_low'])
#    pgrmergelowDict['pointgroup'].append(pg)
#    pgrmergelowDict['Rmerge_I_obs_low'].append(rmergelow)
#    return pgrmergelowDict


def get_crystal_analysis_dict():
    dataDict = {
        'mounted': 0,
        'collected': 0,
        'failed': 0,
        'diffracted': 0,
        'diffracted_2.5': 0
    }
    return dataDict

def get_status_dict():
    statusDict = {}
    return statusDict

def update_crystal_summary(dataDict, field):
    dataDict[field] += 1
    return dataDict


def cell_format_settings(workbook, cif):
    cell_format = workbook.add_format()
    cell_format.set_align('center')
    cell_format.set_align('vcenter')
    cell_format.set_text_wrap()
    if 'FAIL' in cif['status']:
        cell_format.set_font_color('red')
    elif 'OK -' in cif['status']:
        cell_format.set_font_color('orange')
    else:
        cell_format.set_font_color('black')
    return cell_format


def get_unit_cell_string(cif):
    uc = ''
    if 'alpha' in cif:
        uc = (
            '{0!s} '.format(round(float(cif['a']), 1)) +
            '{0!s} '.format(round(float(cif['b']), 1)) +
            '{0!s}\n'.format(round(float(cif['c']), 1)) +
            '{0!s} '.format(round(float(cif['alpha']), 1)) +
            '{0!s} '.format(round(float(cif['beta']), 1)) +
            '{0!s}'.format(round(float(cif['gamma']), 1))
        )
    return uc


def add_row_to_worksheet(workbook, worksheet, sample, cif, row, dozor, cpdID, cpdImg):
    worksheet.set_row(row, 80)
    cell_format = cell_format_settings(workbook, cif)
    unitcell = get_unit_cell_string(cif)
    if cif['pipeline'].startswith('xia2'):
        cif['percent_possible_obs'] = str(round(float(cif['percent_possible_obs'])*100, 1))
    worksheet.write('A' + str(row), sample, cell_format)
    worksheet.write('B' + str(row), cif['sample_appearance'], cell_format)
    worksheet.write('C' + str(row), cif['collection_date'], cell_format)
    worksheet.write('D' + str(row), cif['proposal'], cell_format)
    worksheet.write('E' + str(row), cif['session'], cell_format)
    if dozor:
        worksheet.insert_image('F' + str(row), dozor, {'x_scale': 0.195, 'y_scale': 0.22})
    worksheet.write('G' + str(row), cif['pipeline'], cell_format)
    worksheet.write('H' + str(row), str(cif['reso_low']), cell_format)
    worksheet.write('I' + str(row), str(cif['reso_high']), cell_format)
    worksheet.write('J' + str(row), str(cif['percent_possible_obs']), cell_format)
    worksheet.write('K' + str(row), str(cif['Rmerge_I_obs_low']), cell_format)
    worksheet.write('L' + str(row), str(cif['meanI_over_sigI_obs_high']), cell_format)
    worksheet.write('M' + str(row), str(cif['space_group']), cell_format)
    worksheet.write('N' + str(row), str(unitcell), cell_format)
    worksheet.write('O' + str(row), str(cif['status']), cell_format)
    worksheet.write('P' + str(row), cpdID, cell_format)
    if cpdImg:
        worksheet.insert_image('Q' + str(row), cpdImg, {'x_scale': 0.3, 'y_scale': 0.35})
    worksheet.write('R' + str(row), cif['smiles'], cell_format)
    worksheet.write('S' + str(row), cif['soaktime'], cell_format)
    worksheet.write('T' + str(row), cif['cpd_appearance'], cell_format)

    worksheet.write('U' + str(row), str(cif['initref_software']), cell_format)
    worksheet.write('V' + str(row), str(cif['initref_spacegroup']), cell_format)
    worksheet.write('W' + str(row), '', cell_format)
    worksheet.write('X' + str(row), str(cif['ls_d_res_high']), cell_format)
    worksheet.write('Y' + str(row), str(cif['ls_R_factor_R_work']), cell_format)
    worksheet.write('Z' + str(row), str(cif['ls_R_factor_R_free']), cell_format)
    worksheet.write('AA' + str(row), str(cif['r_bond_refined_d']), cell_format)
    worksheet.write('AB' + str(row), str(cif['r_angle_refined_deg']), cell_format)
    worksheet.write('AC' + str(row), str(cif['blobs']), cell_format)


def get_summary_worksheet(workbook, n_samples):
    worksheet = workbook.add_worksheet('summary')
    merge_format = workbook.add_format({
        'bold': 1,
        'border': 1,
        'align': 'center',
        'valign': 'vcenter',
        'fg_color': 'yellow'})
    merge_format.set_font_size(20)
    worksheet.merge_range('A1:B1', 'Crystal', merge_format)
    worksheet.merge_range('C1:F1', 'Data Collection', merge_format)
    worksheet.merge_range('G1:O1', 'Data Processing', merge_format)
    worksheet.merge_range('P1:T1', 'Compound', merge_format)
    worksheet.merge_range('U1:AC1', 'Initial Refinement', merge_format)
    worksheet.add_table('A2:AC{0!s}'.format(n_samples+2), {'columns': [{'header': 'Sample'},                        # A
                                                                      {'header': 'Appearance'},                     # B
                                                                      {'header': 'Date'},                           # C
                                                                      {'header': 'Proposal'},                       # D
                                                                      {'header': 'Session'},                        # E
                                                                      {'header': 'Dozor'},                          # F
                                                                      {'header': 'ProcPipeline'},                   # G
                                                                      {'header': 'Reso (Low)'},                     # H
                                                                      {'header': 'Reso (High)'},                    # I
                                                                      {'header': 'Completeness\n(Overall)'},        # J
                                                                      {'header': 'Rmerge (Low)'},                   # K
                                                                      {'header': 'I/sig(I) (High)'},                # L
                                                                      {'header': 'ProcSPG'},                        # M
                                                                      {'header': 'Unit Cell'},                      # N
                                                                      {'header': 'Status'},                         # O
                                                                      {'header': 'Compound ID'},                    # P
                                                                      {'header': 'Compound'},                       # Q
                                                                      {'header': 'Smiles'},                         # R
                                                                      {'header': 'Soak time (h)'},                  # S
                                                                      {'header': 'Behaviour'},                      # T
                                                                      {'header': 'RefPipeline'},                    # U
                                                                      {'header': 'RefSPG'},                         # V
                                                                      {'header': 'reference PDB'},                  # W
                                                                      {'header': 'Resolution'},                     # X
                                                                      {'header': 'Rwork'},                          # Y
                                                                      {'header': 'Rfree'},                          # Z
                                                                      {'header': 'rmsd bonds'},                     # AA
                                                                      {'header': 'rmsd angles'},                    # AB
                                                                      {'header': 'blobs'}                           # AC
                                                                      ],
                                                          'header_row': True,
                                                          'autofilter': True})
    worksheet.freeze_panes(2, 1)
    worksheet.set_column('E:E', 17)
    worksheet.set_column('P:P', 12)
    return worksheet


def get_dozor_plot(projectDir, sample, jso, dozor):
    subfolder = jso['proposal'] + '-' + jso['session'] + '-' + jso['run']
    plot = os.path.join(projectDir, '1-process', sample, subfolder, 'images', 'dozor.png')
    if os.path.isfile(plot):
        dozor = plot
    return dozor


def get_compound_image(projectDir, sample, cpdID):
    cpdImg = None
    img = os.path.join(projectDir, '3-compound', sample, cpdID + '.png')
    if os.path.isfile(img):
        cpdImg = img
    return cpdImg


def check_aux_csv_file(auxcsv, sample, cif):
    if os.path.isfile(auxcsv):
        for line in open(auxcsv):
            sampleID = line.split(',')[0]
            statusEntry = line.split(',')[1].replace('\n', '')
            if sampleID == sample:
                cif['status'] = statusEntry
    return cif


def update_status_dict(statusDict, cif):
    status = cif['status']
    if cif['status'] == '':
        status = 'unknown'
    if status not in statusDict:
        statusDict[status] = 0
    statusDict[status] += 1
    return statusDict

def prepare_summary_worksheet_new(logger, workbook, summary_worksheet, projectDir, fragmaxcsv,
                              dataDict, pgDict, pginfoDict, statusDict, auxcsv):
    logger.info('preparing summary worksheet...')
    space_group_dict = {}
    df = pd.read_csv(fragmaxcsv)
    d = df.to_dict(orient='index')
    for m, row in enumerate(d):
        n = m+1
        cif = d[row]

#    for n, line in enumerate(open(fragmaxcsv)):
#        if n == 0:
#            continue
        sample = cif['mounted_crystal_code']
        sample_appearance = cif['crystal_appearance']
        cpdID = cif['compound_code']
        smiles = cif['smiles']
        cpd_appearance = cif['compound_appearance']
        soaktime = "n/a"
        logger.info('current sample: ' + sample)
        dataDict = update_crystal_summary(dataDict, 'mounted')
        ciffile = os.path.join(projectDir, '1-process', sample, 'process.cif')
        mtzfile = os.path.join(projectDir, '1-process', sample, 'process.mtz')
        jsofile = os.path.join(projectDir, '1-process', sample, 'info.json')
        initrefcif = os.path.join(projectDir, '2-initial_refine', sample, 'init.mmcif')
        initrefmtz = os.path.join(projectDir, '2-initial_refine', sample, 'init.mtz')
        initrefpdb = os.path.join(projectDir, '2-initial_refine', sample, 'init.pdb')
        dozor = None
        cpdImg = get_compound_image(projectDir, sample, cpdID)
        if os.path.isfile(jsofile):
            dataDict = update_crystal_summary(dataDict, 'collected')
            jso = get_json_as_dict(jsofile)
            dozor = get_dozor_plot(projectDir, sample, jso, dozor)
        if os.path.isfile(ciffile) and os.path.isfile(mtzfile) and os.path.isfile(jsofile):
#            cif = processlib.cif_info(ciffile)
#            mtz = processlib.mtz_info(mtzfile)
#            cif.update(mtz)
#            cif.update(jso)
#            cif['soaktime'] = soaktime
            dataDict = analyse_resolution(dataDict, cif)
            pgDict = update_point_group_dict(cif, pgDict)
            pginfoDict = update_point_group_info_dict(cif, pginfoDict)
            if os.path.isfile(initrefcif) and os.path.isfile(initrefmtz) and os.path.isfile(initrefpdb):
                cif.update(refinelib.structure_cif_info(initrefcif))
                if cif['initref_spacegroup'] not in space_group_dict:
                    space_group_dict['initref_spacegroup'] = 0
                space_group_dict['initref_spacegroup'] += 1
                blobList = refinelib.find_blobs(initrefmtz, initrefpdb)
                if blobList:
                    cif['blobs'] = str(len(blobList))
                else:
                    cif['blobs'] = '0'
            else:
                cif.update(get_blank_init_refine_cif())
        elif os.path.isfile(jsofile) and not os.path.isfile(ciffile):
#            jso = get_json_as_dict(jsofile)
#            cif = get_semi_blank_cif()
            cif['soaktime'] = soaktime
            cif.update(jso)
        else:
#            cif = get_blank_cif()
            cif['soaktime'] = soaktime
            cif = check_aux_csv_file(auxcsv, sample, cif)
        if sample == 'None':
            cif['status'] = 'FAIL - mounting failed'
        row = n + 2
#        cif['smiles'] = smiles
#        cif['sample_appearance'] = sample_appearance
#        cif['cpd_appearance'] = cpd_appearance
        add_row_to_worksheet(workbook, summary_worksheet, sample, cif, row, dozor, cpdID, cpdImg)
        statusDict = update_status_dict(statusDict, cif)
        #        updated = df['SampleID'] == sample
        updated = df['mounted_crystal_code'] == sample
        uc = str(cif['unitcell']).replace(",", " ")
        cif['unitcell'] = uc
        df_dict = cif
        #        try:
        #            del df_dict['unitcell']
        #        except KeyError:
        #            pass
        #        try:
        #            del df_dict['master']
        #        except KeyError:
        #            pass
        #        print(sample, updated)
        #        print(cif)
        df.loc[updated, df_dict.keys()] = df_dict.values()
    if pginfoDict:
        print_pg_ucv_distribution(pginfoDict)
        print_pg_rmergelow_distribution(pginfoDict)
    if space_group_dict:
        print('space_groups', space_group_dict)
    return statusDict, df


def prepare_summary_worksheet(logger, workbook, summary_worksheet, projectDir, fragmaxcsv,
                              dataDict, pgDict, pginfoDict, statusDict, auxcsv):
    logger.info('preparing summary worksheet...')
    space_group_dict = {}
    df = pd.read_csv(fragmaxcsv)
    df = add_cif_columns_to_df(df)
    for n, line in enumerate(open(fragmaxcsv)):
        if n == 0:
            continue
        sample = line.split(',')[0]
        sample_appearance = line.split(',')[9]
        cpdID = line.split(',')[1]
        smiles = line.split(',')[2]
        cpd_appearance = line.split(',')[10].replace('\n', '').replace('\r', '')
        soaktime = line.split(',')[5]
        logger.info('current sample: ' + sample)
        dataDict = update_crystal_summary(dataDict, 'mounted')
        ciffile = os.path.join(projectDir, '1-process', sample, 'process.cif')
        mtzfile = os.path.join(projectDir, '1-process', sample, 'process.mtz')
        jsofile = os.path.join(projectDir, '1-process', sample, 'info.json')
        initrefcif = os.path.join(projectDir, '2-initial_refine', sample, 'init.mmcif')
        initrefmtz = os.path.join(projectDir, '2-initial_refine', sample, 'init.mtz')
        initrefpdb = os.path.join(projectDir, '2-initial_refine', sample, 'init.pdb')
        dozor = None
        cpdImg = get_compound_image(projectDir, sample, cpdID)
        if os.path.isfile(jsofile):
            dataDict = update_crystal_summary(dataDict, 'collected')
            jso = get_json_as_dict(jsofile)
            dozor = get_dozor_plot(projectDir, sample, jso, dozor)
        if os.path.isfile(ciffile) and os.path.isfile(mtzfile) and os.path.isfile(jsofile):
            cif = processlib.cif_info(ciffile)
            mtz = processlib.mtz_info(mtzfile)
            cif.update(mtz)
            cif.update(jso)
            cif['soaktime'] = soaktime
            dataDict = analyse_resolution(dataDict, cif)
            pgDict = update_point_group_dict(cif, pgDict)
            pginfoDict = update_point_group_info_dict(cif, pginfoDict)
            if os.path.isfile(initrefcif) and os.path.isfile(initrefmtz) and os.path.isfile(initrefpdb):
                cif.update(refinelib.structure_cif_info(initrefcif))
                if cif['initref_spacegroup'] not in space_group_dict:
                    space_group_dict['initref_spacegroup'] = 0
                space_group_dict['initref_spacegroup'] += 1
                blobList = refinelib.find_blobs(initrefmtz, initrefpdb)
                if blobList:
                    cif['blobs'] = str(len(blobList))
                else:
                    cif['blobs'] = '0'
            else:
                cif.update(get_blank_init_refine_cif())
        elif os.path.isfile(jsofile) and not os.path.isfile(ciffile):
            jso = get_json_as_dict(jsofile)
            cif = get_semi_blank_cif()
            cif['soaktime'] = soaktime
            cif.update(jso)
        else:
            cif = get_blank_cif()
            cif['soaktime'] = soaktime
            cif = check_aux_csv_file(auxcsv, sample, cif)
        if sample == 'nan':
            cif['status'] = 'FAIL - mounting failed'
        row = n + 2
        cif['smiles'] = smiles
        cif['sample_appearance'] = sample_appearance
        cif['cpd_appearance'] = cpd_appearance
        add_row_to_worksheet(workbook, summary_worksheet, sample, cif, row, dozor, cpdID, cpdImg)
        statusDict = update_status_dict(statusDict, cif)
#        updated = df['SampleID'] == sample
        updated = df['mounted_crystal_code'] == sample
        uc = cif['unitcell'].replace(",", " ")
        cif['unitcell'] = uc
        df_dict = cif
#        try:
#            del df_dict['unitcell']
#        except KeyError:
#            pass
#        try:
#            del df_dict['master']
#        except KeyError:
#            pass
#        print(sample, updated)
#        print(cif)
        df.loc[updated, df_dict.keys()] = df_dict.values()
    if pginfoDict:
        print_pg_ucv_distribution(pginfoDict)
        print_pg_rmergelow_distribution(pginfoDict)
    if space_group_dict:
        print('space_groups', space_group_dict)
    return statusDict, df


def get_details_worksheet(workbook, n_autoproc_results):
    worksheet = workbook.add_worksheet('details')
    worksheet.add_table('A1:N{0!s}'.format(n_autoproc_results+1), {'columns': [{'header': 'Sample'},
                                                                  {'header': 'Date'},
                                                                  {'header': 'Proposal'},
                                                                  {'header': 'Session'},
                                                                  {'header': 'Run'},
                                                                  {'header': 'Pipeline'},
                                                                  {'header': 'Reso (Low)'},
                                                                  {'header': 'Reso (High)'},
                                                                  {'header': 'Completeness\n(Overall)'},
                                                                  {'header': 'Rmerge (Low)'},
                                                                  {'header': 'I/sig(I) (High)'},
                                                                  {'header': 'Spacegroup'},
                                                                  {'header': 'Unit Cell'},
                                                                  {'header': 'Status'}
                                                                      ]})
    worksheet.freeze_panes(1, 1)
    return worksheet

def get_library_worksheet(workbook, df):
    worksheet = workbook.add_worksheet('libraries')
    lib_df = pd.read_excel("/Users/tobkro/MAXIV/FragMAX/LP3/FragMAX_crystal_preparation/database/lookup_tables.xlsx", 6)
    ld = {}
    x = 1
    y = 1
    for index, row in lib_df.iterrows():
        cpd_id = row['compound_code']
        lib_name = row['library_name']
        if lib_name not in ld:
            ld[lib_name] = []
            y += 1
            x = 1
        if cpd_id not in ld[lib_name]:
            ld[lib_name].append(cpd_id)
            if x == 27:
                x = 1
                y += 1
            col = number_to_column(str(x))
            dt = df.loc[df['compound_code'] == cpd_id]
#            dt = df.loc[df['CompoundID'] == cpd_id]
            for i, r in dt.iterrows():
                try:
                    if float(r['ls_d_res_high']) < 2.5:
                        cell_format = workbook.add_format({'bg_color': '00FF00'})
                    else:
                        cell_format = workbook.add_format({'bg_color': '006666'})
                except ValueError:
#                    cell_format = workbook.add_format({'bg_color': 'red'})
#                    print('>>>', r['mounted_crystal_code'], r['ls_d_res_high'], r['sample_appearance'])
                    if 'bad' in r['sample_appearance'].lower():
                        cell_format = workbook.add_format({'bg_color': 'yellow'})
                    elif r['mounted_crystal_code'] == 'None':
                        cell_format = workbook.add_format({'bg_color': 'magenta'})  # something went wrong during mount
                    elif r['mounted_crystal_code'] != 'None' and r['status'] == "FAIL - no processing result":
                        cell_format = workbook.add_format({'bg_color': 'red'})
                    elif r['mounted_crystal_code'] != 'None' and r['status'].replace(' ', '') == '':    # prepared but not collected
                        cell_format = workbook.add_format({'bg_color': 'cyan'})
#                    elif r['status'].replace(' ', '') == '':
#                        cell_format = workbook.add_format({'bg_color': 'cyan'})
                    else:
                        cell_format = workbook.add_format({'bg_color': 'brown'})    # all other failure
            if dt.empty:
                cell_format = workbook.add_format({'bg_color': 'white'})

            worksheet.write('{0!s}{1!s}'.format(col, str(y)), cpd_id, cell_format)
            x += 1
        else:
            continue
    return worksheet


def prepare_details_worksheet(logger, workbook, details_worksheet, projectDir, fragmaxcsv):
    logger.info('preparing details worksheet...')
    row = 1 # since it has a header
    color_one = '#f2ebeb'
    color_two = '#d9d2d2'
    color = color_one
    for n, line in enumerate(open(fragmaxcsv)):
        sample = line.split(',')[0]
        logger.info('current sample: ' + sample)
        for cif in glob.glob(os.path.join(projectDir, '1-process', sample, '*', '*', 'process.cif')):
            ciffile = cif
            mtzfile = cif.replace('process.cif', 'process.mtz')
            jsofile = cif.replace('process.cif', 'info.json')
            if os.path.isfile(ciffile) and os.path.isfile(mtzfile) and os.path.isfile(jsofile):
                cif = processlib.cif_info(ciffile)
                mtz = processlib.mtz_info(mtzfile)
                jso = get_json_as_dict(jsofile)
                cif.update(mtz)
                cif.update(jso)
                row += 1
                add_row_to_details_worksheet(workbook, details_worksheet, sample, cif, row, color)
        if (n % 2) == 0:
            color = color_one
        else:
            color = color_two


def cell_format_details_settings(workbook, color):
    cell_format = workbook.add_format()
    cell_format.set_bg_color(color)
    return cell_format


def add_row_to_details_worksheet(workbook, worksheet, sample, cif, row, color):
    cell_format = cell_format_details_settings(workbook, color)
#    if cif['pipeline'].startswith('xia2'):
#        cif['percent_possible_obs'] = str(round(float(cif['percent_possible_obs'])*100, 1))
    worksheet.write('A' + str(row), sample, cell_format)
    worksheet.write('B' + str(row), cif['collection_date'], cell_format)
    worksheet.write('C' + str(row), cif['proposal'], cell_format)
    worksheet.write('D' + str(row), cif['session'], cell_format)
    worksheet.write('E' + str(row), cif['run'], cell_format)
    worksheet.write('F' + str(row), cif['pipeline'], cell_format)
    worksheet.write('G' + str(row), cif['reso_low'], cell_format)
    worksheet.write('H' + str(row), cif['reso_high'], cell_format)
    worksheet.write('I' + str(row), cif['percent_possible_obs'], cell_format)
    worksheet.write('J' + str(row), cif['Rmerge_I_obs_low'], cell_format)
    worksheet.write('K' + str(row), cif['meanI_over_sigI_obs_high'], cell_format)
    worksheet.write('L' + str(row), cif['space_group'], cell_format)
    worksheet.write('M' + str(row), cif['unitcell'], cell_format)
    worksheet.write('N' + str(row), cif['status'], cell_format)


def get_pg_ucv_worksheet(workbook):
    worksheet = workbook.add_worksheet('PG-UCV')
    return worksheet


def prepare_get_pg_ucv_worksheet(pg_ucv_worksheet):
    if os.path.isfile('pg_ucv_distribution.png'):
        pg_ucv_worksheet.insert_image('A1', 'pg_ucv_distribution.png')


def print_pg_ucv_distribution(pgucvDict):
    df = pd.DataFrame(pgucvDict)
    pg_plot = sns.catplot(x="pointgroup", y="unitcell_volume", data=df)
#    pg_plot = sns.catplot(x="pointgroup", y="unitcell_volume", data=df, kind="violin")
    plt.savefig('pg_ucv_distribution.png', dpi=300)


def get_pg_rmergelow_worksheet(workbook):
    worksheet = workbook.add_worksheet('PG-RmergeLow')
    return worksheet


def prepare_get_pg_rmergelow_worksheet(pg_ucv_worksheet):
    if os.path.isfile('pg_ucv_distribution.png'):
        pg_ucv_worksheet.insert_image('A1', 'pg_rmergelow_distribution.png')


def print_pg_rmergelow_distribution(pgucvDict):
    df = pd.DataFrame(pgucvDict)
    pg_plot = sns.catplot(x="pointgroup", y="Rmerge_I_obs_low", data=df)
#    pg_plot = sns.catplot(x="pointgroup", y="Rmerge_I_obs_low", data=df, kind="violin")
    plt.savefig('pg_rmergelow_distribution.png', dpi=300)


def analyse_resolution(dataDict, cif):
    if float(cif['reso_high']) < 2.5:
        dataDict = update_crystal_summary(dataDict, 'diffracted_2.5')
        dataDict = update_crystal_summary(dataDict, 'diffracted')
    else:
        dataDict = update_crystal_summary(dataDict, 'diffracted')
    return dataDict


def get_chart_data_sheet(workbook):
    chart_data = workbook.add_worksheet('chart data')
    return chart_data


def prepare_crystal_chart(workbook, chart_data, dataDict):
    n = 0
    for n, field in enumerate(dataDict):
        chart_data.write(number_to_column(n + 1) + '1', field)
        chart_data.write(number_to_column(n + 1) + '2', dataDict[field])
    chartsheet = workbook.add_chartsheet('crystals')
    chart = workbook.add_chart({'type': 'column'})

    chart.add_series({
        'categories': ['chart data', 0, 0, 0, n],
        'values': ['chart data', 1, 0, 1, n],
        'line': {'color': 'black'},
    })
    chart.set_title({'name': 'Crystals'})
    chart.set_legend({'none': True})
    chartsheet.set_chart(chart)
#    chartsheet.set_zoom(200)


def prepare_laue_group_chart(workbook, chart_data, pgDict):
    n = 0
    for n, field in enumerate(pgDict):
        chart_data.write(number_to_column(n + 1) + '4', field)
        chart_data.write(number_to_column(n + 1) + '5', pgDict[field])
    chartsheet_pg = workbook.add_chartsheet('pointgroup')
    chart_pg = workbook.add_chart({'type': 'pie'})

    chart_pg.add_series({
        'categories': ['chart data', 3, 0, 3, n],
        'values': ['chart data', 4, 0, 4, n],
        'line': {'color': 'black'},
        'data_labels': {'value': True},
    })
    chart_pg.set_title({'name': 'point group'})
    chart_pg.set_legend({'font': {'size': 20, 'bold': True}})
    #    chart_pg.set_legend({'none': True})
    chartsheet_pg.set_chart(chart_pg)


def prepare_status_chart(workbook, chart_data, statusDict):
    n = 0
    for n, field in enumerate(statusDict):
        chart_data.write(number_to_column(n + 1) + '7', field)
        chart_data.write(number_to_column(n + 1) + '8', statusDict[field])
    chartsheet_status = workbook.add_chartsheet('status')
    chart_status = workbook.add_chart({'type': 'pie'})

    chart_status.add_series({
        'categories': ['chart data', 6, 0, 6, n],
        'values': ['chart data', 7, 0, 7, n],
        'line': {'color': 'black'},
        'data_labels': {'value': True},
    })
    chart_status.set_title({'name': 'status'})
    chart_status.set_legend({'font': {'size': 20, 'bold': True}})
    #    chart_pg.set_legend({'none': True})
    chartsheet_status.set_chart(chart_status)
