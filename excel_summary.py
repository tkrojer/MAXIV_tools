# Copyright (c) 2022, Tobias Krojer, MAX IV Laboratory
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import sys
import os
import xlsxwriter
from PIL import Image
import getopt
import glob


sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'lib')))
import processlib
import summarylib

#with Image.open("hsDHS-x0001_1_1.snapshot.jpeg") as img:
#    width_100 = img.width
#    height_100 = img.height
#    print(width_100, height_100)
##    print(img.info['dpi'])

# save a smaller version of the image
#width_30 = int(round(width_100 * 0.3, 0))
#img = Image.open('local_100_perc.png')
#wpercent = (width_30/float(width_100))
#hsize = int((float(height_100)*float(wpercent)))
#img = img.resize((width_30,hsize), Image.ANTIALIAS)
#img.save('local_30_perc.png')


#workbook = xlsxwriter.Workbook('image.xlsx')
#worksheet = workbook.add_worksheet()

# worksheets
# Data collection - selected
# Data collection - all






#worksheet.add_table('A1:C3', {'columns': [{'header': 'Sample'},
 #                                         {'header': 'Dozor'},
 #                                         {'header': 'tbc'}
 #                                         ]})

#wrap_format = workbook.add_format({'text_wrap': True})

# locks cells; see https://xlsxwriter.readthedocs.io/example_protection.html
#worksheet.protect()

#worksheet.write('A1', 'sample 1')
#worksheet.write('A2', 'sample 1', wrap_format)
#worksheet.set_column('B:B', 80)
#worksheet.set_row(0, 80)
#worksheet.insert_image('B2', 'hsDHS-x0001_1_1.snapshot.jpeg', {'object_position': 1})
#worksheet.insert_image('C1', 'hsDHS-x0001_1_1.snapshot.jpeg')
#worksheet.insert_image('B8', 'logo.png')

#workbook.close()


# input
# projectDir
# fragmaxcsv
# auxcsv - commma separated file with additinal info like "Fail - broken loop", "Fail - empty loop"

#def read_project_information():


def add_row_to_worksheet(worksheet, cell_format, sample, cif, row):
    worksheet.write('A' + str(row), sample, cell_format)
    worksheet.write('B' + str(row), cif['collection_date'], cell_format)
    worksheet.write('C' + str(row), cif['proposal'], cell_format)
    worksheet.write('D' + str(row), cif['session'], cell_format)
    worksheet.write('E' + str(row), cif['reso_low'], cell_format)
    worksheet.write('F' + str(row), cif['reso_high'], cell_format)
    worksheet.write('G' + str(row), cif['status'], cell_format)


def get_crystal_analysis_dict():
    dataDict = {
        'mounted': 0,
        'collected': 0,
        'failed': 0,
        'diffracted': 0,
        'diffracted_2.5': 0
    }
    return dataDict


def get_point_group_dict():
    pgDict = {}
    return pgDict


def update_point_group_dict(cif, pgDict):
    pg = cif['lattice'] + cif['point_group']
    if pg not in pgDict:
        pgDict[pg] = 0
    pgDict[pg] += 1
    return pgDict



def update_crystal_summary(dataDict, field):
    dataDict[field] += 1
    return dataDict

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
        '11': 'K'
    }
    return columnDict[str(n)]

def parse_project_directory(projectDir, fragmaxcsv):
    dataDict = get_crystal_analysis_dict()
    pgDict = get_point_group_dict()
    workbook = summarylib.init_workbook()
    n_samples = summarylib.get_n_samples(fragmaxcsv)
    worksheet = workbook.add_worksheet('summary')
    cell_format = workbook.add_format()
    cell_format.set_align('center')
    cell_format.set_align('vcenter')
    cell_format.set_text_wrap()

    merge_format = workbook.add_format({
        'bold': 1,
        'border': 1,
        'align': 'center',
        'valign': 'vcenter',
        'fg_color': 'yellow'})
    merge_format.set_font_size(20)

    worksheet.merge_range('A1:D1', 'Data Collection', merge_format)
    worksheet.merge_range('E1:K1', 'Data Processing', merge_format)

    worksheet.add_table('A2:G{0!s}'.format(n_samples+2), {'columns': [{'header': 'Sample'},
                                                                      {'header': 'Date'},
                                                                      {'header': 'Proposal'},
                                                                      {'header': 'Session'},
                                                                      {'header': 'Pipeline'},
                                                                      {'header': 'Reso (Low)'},
                                                                      {'header': 'Reso (High)'},
                                                                      {'header': 'Completeness\n(Overall)'},
                                                                      {'header': 'Rmerge (Low)'},
                                                                      {'header': 'I/sig(I) (High)'},
                                                                      {'header': 'Status'}
                                                                      ]})

    worksheet.freeze_panes(2, 1)



    for n, line in enumerate(open(fragmaxcsv)):
        sample = line.split(',')[0]
        dataDict = update_crystal_summary(dataDict, 'mounted')
        ciffile = os.path.join(projectDir, '1-process', sample, 'process.cif')
        mtzfile = os.path.join(projectDir, '1-process', sample, 'process.mtz')
        jsofile = os.path.join(projectDir, '1-process', sample, 'info.json')
        if os.path.isfile(jsofile):
            dataDict = update_crystal_summary(dataDict, 'collected')
        if os.path.isfile(ciffile) and os.path.isfile(mtzfile) and os.path.isfile(jsofile):
            cif = processlib.cif_info(ciffile)
            mtz = processlib.mtz_info(mtzfile)
            jso = summarylib.get_json_as_dict(jsofile)
            cif.update(mtz)
            cif.update(jso)
            dataDict = analyse_resolution(dataDict, cif)
            pgDict = update_point_group_dict(cif, pgDict)
        elif os.path.isfile(jsofile) and not os.path.isfile(ciffile):
            jso = summarylib.get_json_as_dict(jsofile)
            cif = summarylib.get_semi_blank_cif()
            cif.update(jso)
        else:
            cif = summarylib.get_blank_cif()
        row = n + 3
        add_row_to_worksheet(worksheet, cell_format, sample, cif, row)

#    print(dataDict)

    chart_data = workbook.add_worksheet('chart data')
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
    chartsheet.set_zoom(200)


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
    })
    print(n+1)
    chart_pg.set_title({'name': 'point group'})
    chart_pg.set_legend({'font': {'size': 20, 'bold': True}})
#    chart_pg.set_legend({'none': True})
    chartsheet_pg.set_chart(chart_pg)



    print(pgDict)



    workbook.close()


def analyse_resolution(dataDict, cif):
    if float(cif['reso_high']) < 2.5:
        dataDict = update_crystal_summary(dataDict, 'diffracted_2.5')
        dataDict = update_crystal_summary(dataDict, 'diffracted')
    else:
        dataDict = update_crystal_summary(dataDict, 'diffracted')
    return dataDict

def main(argv):
    projectDir = ''
    fragmaxcsv = ''
    auxcsv = ''
    logger = processlib.init_logger('excel_summary.log')
    processlib.start_logging(logger, 'excel_summary.py')

    try:
        opts, args = getopt.getopt(argv,"p:f:a:h",["project=", "fragmax=", "auxcsvn=", "help"])
    except getopt.GetoptError:
        processlib.usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            processlib.usage()
            sys.exit(2)
        elif opt in ("-p", "--project"):
            projectDir = os.path.abspath(arg)
        elif opt in ("-f", "--fragmaxcsv"):
            fragmaxcsv = os.path.abspath(arg)
        elif opt in ("-a", "--auxcsv"):
            auxcsv = os.path.abspath(arg)

#    processlib.report_parameters(logger, processDir, projectDir, fragmaxcsv, select, select_criterion, overwrite)
#    checks_passed = processlib.run_checks(logger, processDir, projectDir, fragmaxcsv, select, select_criterion)

    parse_project_directory(projectDir, fragmaxcsv)

#    if checks_passed:
#        processlib.check_if_to_continue(logger)
#        if select:
#            processlib.start_select_results(logger)
#            select_results(logger, projectDir, select_criterion, overwrite)
#        else:
#            processlib.start_get_autoprocessing_results(logger)
#            get_autoprocessing_results(logger, processDir, projectDir, fragmaxcsv)
#    else:
#        logger.error('cannot continue; check error messages above and use -h option to get more information')
#        logger.info('===================================================================================')

if __name__ == '__main__':
    main(sys.argv[1:])
