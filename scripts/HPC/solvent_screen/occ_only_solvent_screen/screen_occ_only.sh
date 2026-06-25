#!/bin/bash -l
#SBATCH -p NAME_QUEUE
#SBATCH -c NUM_CORES
#SBATCH --array=1-179

# Use a default input deck unless specified on the command-line
if [ -n "$1" ]; then
  CIF=$1
else
  # Need to be in the current dir
    CIF=MY_CIF.cif
fi

# Activate python environment if using occpy
source PYTHON_ENVIRONMENT_PATH

# Extract Solvent based on JobArrayID
TASK_ID=${SLURM_ARRAY_TASK_ID}
SOLVENT_FILE="solvents.txt"
SOLVENT=$(sed -n "${TASK_ID}p" ${SOLVENT_FILE} | awk '{$1=$1; print}')

SIM_DIR=$(pwd)
cif_base=$(basename "$CIF" .cif)

# Debug: Print TASK_ID and SOLVENT
# echo "TASK_ID: ${TASK_ID}"
# echo "SOLVENT: ${SOLVENT}"
# echo "SIM_DIR: ${SIM_DIR}"

# Run OCC with specific solvent
echo "Computing solvation ${CIF} in ${SOLVENT}"
occpy cg ${CIF} --radius=30 --cg-radius=3.8 --threads=${SLURM_CPUS_PER_TASK} --solvent="${SOLVENT}" > "${CIF%.cif}.${SOLVENT}.stdout"
echo "Done with OCC solvation computation"


