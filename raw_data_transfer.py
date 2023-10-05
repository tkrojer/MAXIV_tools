import os
import paramiko
import json
import getpass
import sys
import getopt

def local_sample_folder_exists(local_dir, sample):
    sample_folder_exists = False
    os.chdir(local_dir)
    if not os.path.isdir(sample):
        print('creating sample directory {0!s} in {1!s}'.format(sample, local_dir))
        os.mkdir(sample)
    if os.path.isdir(os.path.join(local_dir, sample)):
        print('SUCCESS: local sample directory exists')
        sample_folder_exists = True
    else:
        print('ERROR: local sample directory does not exist')
    return sample_folder_exists


def copy_diffraction_data(username, password, sftp_server, project_dir, local_dir, sample_list):
    print('connecting to {0!s}...'.format(sftp_server))
    ssh_client = paramiko.SSHClient()
    ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh_client.connect(sftp_server, username=username, password=password)
    sftp_client=ssh_client.open_sftp()
    for sample in sample_list:
        print('current sample: {0!s}'.format(sample))
        if local_sample_folder_exists(local_dir, sample):
            os.chdir(os.path.join(local_dir, sample))
            try:
                remote_file = sftp_client.open(os.path.join(project_dir, '1-process', sample, 'info.json'))
                data = json.load(remote_file)
                master_file = data['master']
                remote_path = master_file[:master_file.rfind('/')]
                run = master_file.replace(remote_path, '').replace('master.h5', '').replace('/', '')
                for files in sftp_client.listdir(remote_path):
                    file_name = files[files.rfind('/')+1:]
                    if run in file_name:
                        data_file = os.path.join(remote_path, file_name)
                        print('copying {0!s}...'.format(file_name))
                        sftp_client.get(data_file, os.path.join(local_dir, sample, file_name))
            except FileNotFoundError:
                print('ERROR: info.json does not exist in {0!s}'.format(os.path.join(project_dir, '1-process', sample)))

    ssh_client.close()
    del ssh_client, sftp_client
    exit_code = 0
    return exit_code


def get_sample_list(dataset_file):
    sample_list = []
    for line in open(dataset_file):
        sample_list.append(line.replace(' ','').replace('\n','').replace('\r', ''))
    return sample_list


def main(argv):
    sftp_server = 'sftp.maxiv.lu.se'
    project_dir = None
    local_dir = None
    dataset_file = None

    try:
        opts, args = getopt.getopt(argv, "p:d:l:h",["project_dir=", "dataset_file=", "local_dir=", "help"])
    except getopt.GetoptError:
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print('ccp4-python raw_data_transfer.py'
                  ' -p <fragmax_project_directory>'
                  ' -d <dataset_text_file>'
                  ' -l <local_ddirectory> ')
            sys.exit(2)
        elif opt in ("-p", "--project_dir"):
            project_dir = arg
        elif opt in ("-l", "--local_dir"):
            local_dir = os.path.abspath(arg)
        elif opt in ("-d", "--dataset_file"):
            dataset_file = arg

    if project_dir and local_dir and os.path.isfile(dataset_file):
        username = input('username: ')
        password = getpass.getpass('password: ')
        sample_list = get_sample_list(dataset_file)
        exit_code = copy_diffraction_data(username, password, sftp_server, project_dir, local_dir, sample_list)
        sys.exit(exit_code)
    else:
        print('INPUT ERROR: run "ccp4-python raw_data_transfer.py -h"')

if __name__ == '__main__':
    main(sys.argv[1:])
