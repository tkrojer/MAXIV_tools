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
        'space_group': ''
    }
    return cif


def get_semi_blank_cif():
    cif = {
        'reso_low': '',
        'reso_high': '',
        'pipeline': '',
        'percent_possible_obs': '',
        'Rmerge_I_obs_low': '',
        'meanI_over_sigI_obs_high': '',
        'unitcell': '',
        'space_group': ''
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
    pg = cif['lattice'] + cif['point_group']
    if pg not in pgDict:
        pgDict[pg] = 0
    pgDict[pg] += 1
    return pgDict


def get_point_group_ucv_dict():
    pgucvDict = {
        'pointgroup': [],
        'unitcell_volume': []
    }
    return pgucvDict


def update_point_group_ucv_dict(cif, pgucvDict):
    pg = cif['lattice'] + cif['point_group']
    ucv = int(cif['unitcell_volume'])
    pgucvDict['pointgroup'].append(pg)
    pgucvDict['unitcell_volume'].append(ucv)
    return pgucvDict


def get_crystal_analysis_dict():
    dataDict = {
        'mounted': 0,
        'collected': 0,
        'failed': 0,
        'diffracted': 0,
        'diffracted_2.5': 0
    }
    return dataDict


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
    worksheet.write('B' + str(row), cif['collection_date'], cell_format)
    worksheet.write('C' + str(row), cif['proposal'], cell_format)
    worksheet.write('D' + str(row), cif['session'], cell_format)
    if dozor:
        worksheet.insert_image('E' + str(row), dozor, {'x_scale': 0.195, 'y_scale': 0.22})
    worksheet.write('F' + str(row), cif['pipeline'], cell_format)
    worksheet.write('G' + str(row), cif['reso_low'], cell_format)
    worksheet.write('H' + str(row), cif['reso_high'], cell_format)
    worksheet.write('I' + str(row), cif['percent_possible_obs'], cell_format)
    worksheet.write('J' + str(row), cif['Rmerge_I_obs_low'], cell_format)
    worksheet.write('K' + str(row), cif['meanI_over_sigI_obs_high'], cell_format)
    worksheet.write('L' + str(row), cif['space_group'], cell_format)
    worksheet.write('M' + str(row), unitcell, cell_format)
    worksheet.write('N' + str(row), cif['status'], cell_format)
    worksheet.write('O' + str(row), cpdID, cell_format)
    if cpdImg:
        worksheet.insert_image('P' + str(row), cpdImg, {'x_scale': 0.3, 'y_scale': 0.35})


def get_summary_worksheet(workbook, n_samples):
    worksheet = workbook.add_worksheet('summary')
    merge_format = workbook.add_format({
        'bold': 1,
        'border': 1,
        'align': 'center',
        'valign': 'vcenter',
        'fg_color': 'yellow'})
    merge_format.set_font_size(20)
    worksheet.merge_range('A1:E1', 'Data Collection', merge_format)
    worksheet.merge_range('F1:N1', 'Data Processing', merge_format)
    worksheet.merge_range('O1:P1', 'Compound', merge_format)
    worksheet.add_table('A2:P{0!s}'.format(n_samples+2), {'columns': [{'header': 'Sample'},
                                                                      {'header': 'Date'},
                                                                      {'header': 'Proposal'},
                                                                      {'header': 'Session'},
                                                                      {'header': 'Dozor'},
                                                                      {'header': 'Pipeline'},
                                                                      {'header': 'Reso (Low)'},
                                                                      {'header': 'Reso (High)'},
                                                                      {'header': 'Completeness\n(Overall)'},
                                                                      {'header': 'Rmerge (Low)'},
                                                                      {'header': 'I/sig(I) (High)'},
                                                                      {'header': 'Spacegroup'},
                                                                      {'header': 'Unit Cell'},
                                                                      {'header': 'Status'},
                                                                      {'header': 'Compound ID'},
                                                                      {'header': 'Compound'}
                                                                      ]})
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


def prepare_summary_worksheet(workbook, summary_worksheet, projectDir, fragmaxcsv, dataDict, pgDict, pgucvDict):
    for n, line in enumerate(open(fragmaxcsv)):
        sample = line.split(',')[0]
        cpdID = line.split(',')[1]
        dataDict = update_crystal_summary(dataDict, 'mounted')
        ciffile = os.path.join(projectDir, '1-process', sample, 'process.cif')
        mtzfile = os.path.join(projectDir, '1-process', sample, 'process.mtz')
        jsofile = os.path.join(projectDir, '1-process', sample, 'info.json')
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
            dataDict = analyse_resolution(dataDict, cif)
            pgDict = update_point_group_dict(cif, pgDict)
            pgucvDict = update_point_group_ucv_dict(cif, pgucvDict)
        elif os.path.isfile(jsofile) and not os.path.isfile(ciffile):
            jso = get_json_as_dict(jsofile)
            cif = get_semi_blank_cif()
            cif.update(jso)
        else:
            cif = get_blank_cif()
        row = n + 3
        add_row_to_worksheet(workbook, summary_worksheet, sample, cif, row, dozor, cpdID, cpdImg)
    if pgucvDict:
        print_pg_ucv_distribution(pgucvDict)


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


def prepare_details_worksheet(workbook, details_worksheet, projectDir, fragmaxcsv):
    row = 1 # since it has a header
    color_one = '#f2ebeb'
    color_two = '#d9d2d2'
    color = color_one
    for n, line in enumerate(open(fragmaxcsv)):
        sample = line.split(',')[0]
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
    plt.savefig('pg_ucv_distribution.png', dpi=300)


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


#def prepare_dozor_thumbnail()