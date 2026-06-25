#!/bin/bash -l
#SBATCH -p QUEUE_NAME
#SBATCH -c NUM_CORES

source PATH_TO_PYTHON_ENV

# Use a default input deck unless specified on the command-line
if [ -n "$1" ]; then
  CIF=$1
else
  # Need to be in the current dir
    CIF=MY_CIF.cif
fi

occpy cg ${CIF} --radius=30 --cg-radius=3.8 --threads=${SLURM_CPUS_PER_TASK} > ${CIF%.cif}.stdout

