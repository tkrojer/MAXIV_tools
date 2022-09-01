import xlsxwriter
from datetime import datetime
import json


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
        'meanI_over_sigI_obs_high': ''    }
    return cif


def get_semi_blank_cif():
    cif = {
        'reso_low': '',
        'reso_high': '',
        'pipeline': '',
        'percent_possible_obs': '',
        'Rmerge_I_obs_low': '',
        'meanI_over_sigI_obs_high': ''
    }
    return cif


def get_json_as_dict(jsonfile):
    with open(jsonfile) as d:
        dictData = json.load(d)
    return dictData