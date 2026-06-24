#!/bin/bash -l
#$ -cwd
#$ -pe smp.pe 12

# Use a default input deck unless specified on the command-line
if [ -n "$1" ]; then
  CIF=$1
else
  # Need to be in the current dir
    CIF=MY_CIF.cif
fi

export OCC_DATA_PATH="YOUR_PATH/occ/share/"
occ=YOUR_PATH/occ/bin/occ-0.6.12-linux-x86_64-static/bin/occ
${occ} cg ${CIF} --radius=30 --cg-radius=3.8 --threads=${NSLOTS} > ${CIF%.cif}.stdout

cp ${CIF%.cif}_cg.txt STRUCTURE_FILE.txt
