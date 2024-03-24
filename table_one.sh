#!/bin/bash

module load gopresto CCP4

# Initialize variables for arguments
model=""
cif=""
xml=""

# Parse arguments
while getopts "m:c:x:" opt; do
  case $opt in
    m) model=$OPTARG ;;
    s) cif=$OPTARG ;;
    o) xml=$OPTARG ;;
    \?) echo "Invalid option -$OPTARG" >&2 ;;
  esac
done

# Call the Python script with the parsed arguments
ccp4-python /data/staff/biomax/tobias/software/MAXIV_tools/table_one.py -m "$model" -c "$cif" -x "$xml"
