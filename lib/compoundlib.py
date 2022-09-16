import glob
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from datetime import datetime


def create_sample_folder(logger, projectDir, sample):
    os.chdir(os.path.join(projectDir, '3-compound'))
    if os.path.isdir(sample):
        logger.warning('sample folder exists')
    else:
        logger.info('creating folder {0!s} in 3-compound')
        os.mkdir(sample)


def empty_sample_folder(logger, projectDir, sample):
    os.chdir(os.path.join(projectDir, '3-compound', sample))
    logger.warning('removing all contents from folder because "--overwrite" was selected')
    os.system('rm -fr *')


def make_compound_png(logger, projectDir, sample, cpdID, smiles):
    os.chdir(os.path.join(projectDir, '3-compound', sample))
    if os.path.isfile("{0!s}.png".format(cpdID)):
        logger.warning("{0!s}.png exisits; skipping...".format(cpdID))
    else:
        logger.info("creating {0!s}.png with rdkit".format(cpdID))
        mol = Chem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(mol)
        Draw.MolToFile(mol, "{0!s}.png".format(cpdID))
        if os.path.isfile("{0!s}.png".format(cpdID)):
            logger.info("{0!s}.png successfully created".format(cpdID))
        else:
            logger.error("creation of {0!s}.png failed".format(cpdID))


def backup_existing_script_files(script):
    if os.path.isfile(script):
        now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        os.system('/bin/mv {0!s} {1!s}_{2!s}'.format(script, script, now))


def ligand_files_exist(logger, cpdID):
    files_exist = False
    if os.path.isfile(cpdID + '.cif') or os.path.isfile(cpdID + '.cif'):
        logger.warning('CIF and PDB files exist for compound; skipping...')
        files_exist = True
    return files_exist


def modules_to_load(restraints_program):
    if restraints_program == 'acedrg':
        module = 'module load gopresto CCP4\n'
    elif restraints_program == 'grade':
        module = 'module load gopresto BUSTER\n'
    elif restraints_program == 'elbow':
        module = 'module load gopresto Phenix\n'
    return module


def maxiv_header(restraints_program):
    header = (
        '#!/bin/bash\n'
        '#SBATCH --time=10:00:00\n'
        '#SBATCH --job-name={0!s}\n'.format(restraints_program) +
        '#SBATCH --cpus-per-task=1\n'
    )
    return header


def restraints_program_cmd(restraints_program, ligandID, smiles):
    if restraints_program == 'acedrg':
        cmd = 'acedrg --res LIG -i "{0!s}" -o {1!s}'.format(smiles, ligandID)
    elif restraints_program == 'grade':
        cmd = ''
    elif restraints_program == 'elbow':
        cmd = ''
    return cmd


def prepare_script_for_maxiv(restraints_program, projectDir, sampleID, subdirectory, ligandID, smiles):
    os.chdir(os.path.join(projectDir, sampleID, subdirectory))
    cmd = maxiv_header()
    cmd += maxiv_header(restraints_program)
    cmd += 'cd {0!s}\n'.format(os.path.join(projectDir, sampleID, subdirectory))
    cmd += restraints_program_cmd(restraints_program, ligandID, smiles)
    f = open('{0!s}.sh'.format(restraints_program), 'w')
    f.write(cmd)
    f.close()


def prepare_script_for_restrains_generation(logger, projectDir, sample, cpdID, smiles, software, submitList):
    os.chdir(os.path.join(projectDir, 'tmp'))
    if not ligand_files_exist(logger, cpdID):
        logger.info('creating script for {0!s}'.format(software))
        script = software + '_' + sample + '.sh'
        backup_existing_script_files(script)
        cmd = maxiv_header(software)
        cmd += modules_to_load(software)
        cmd += 'cd {0!s}\n'.format(os.path.join(projectDir, '3-compound', sample))
        cmd += restraints_program_cmd(software, cpdID, smiles)
        f = open(script, 'w')
        f.write(cmd)
        f.close()
        submitList.append(script)
    return submitList


def submit_jobs_to_cluster(logger, projectDir, submitList):
    logger.info('submitting jobs to MAXIV cluster...')
    os.chdir(os.path.join(projectDir, 'tmp'))
    for script in submitList:
        logger.info('submitting ' + script)
        os.system('sbatch ' + script)


def usage():
    usage = (
        '\n'
        'usage:\n'
        'ccp4-python 3-compound.py -p <project_dir> -f <fragmax_csv_file>\n'
        '\n'
        'additional command line options:\n'
        '--project, -p\n'
        '    project directory\n'
        '--fragmax, -f\n'
        '    fragmax summary csv file\n'
        '--overwrite, -o\n'
        '    flag to overwrite selected files\n'
    )
    print(usage)



