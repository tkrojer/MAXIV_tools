#!/bin/bash

module load gopresto CCP4

# Initialize variables for arguments
model=""
sf=""
output=""

# Parse arguments
while getopts "i:o:x:" opt; do
  case $opt in
    m) model=$OPTARG ;;
    s) sf=$OPTARG ;;
    o) output=$OPTARG ;;
    \?) echo "Invalid option -$OPTARG" >&2 ;;
  esac
done

# Call the Python script with the parsed arguments
ccp4-python /data/staff/biomax/tobias/software/MAXIV_tools/pdb_validate_cif.py -m "$model" -s "$sf" -o "$output"
