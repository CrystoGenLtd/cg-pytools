#!/bin/bash -l
#$ -cwd
#$ -pe smp.pe 4
#$ -t 1-179

# Use a default input deck unless specified on the command-line
if [ -n "$1" ]; then
  CIF=$1
else
  # Need to be in the current dir
    CIF=MY_CIF.cif
fi

# Set environment variables, path to CIF, and OCC/CG executables  
CG_KEY="YOUR_PATH/crystogen-1.2.0-linux-x86_64/bin/CrystoGen.key"
CG_EXE="YOUR_PATH/crystogen-1.2.0-linux-x86_64/bin/crystogen"

# Extract Solvent based on JobArrayID
TASK_ID=${SGE_TASK_ID}
SOLVENT_FILE="solvents.txt"
SOLVENT=$(sed -n "${TASK_ID}p" ${SOLVENT_FILE} | awk '{$1=$1; print}')

SIM_DIR=$(pwd)
cif_base=$(basename "$CIF" .cif)
simulation_path=$(realpath "$SIM_DIR")

if [ $TASK_ID -eq 1 ]
then
# Copy the facets all file into parent folder and call it COLOUR.txt
# add number of facets to first line to make it work
    echo "MAKING TEMPLATE COLOUR FILE"
    colour=$(find . -type f -name "*all.txt" -not -empty -print -quit)
    lines=$(wc -l < $colour)
    rm COLOUR.txt
    echo "  ${lines}" >> COLOUR.txt
    cat $colour >> COLOUR.txt
fi

# Path to the directory where the simulation should run
solvent="${SOLVENT// /+}"
simulation_path="$SIM_DIR"

# Create solvent-specific folder
solvent_folder="$simulation_path/solvent_$solvent"

# Modify and copy input.txt file to solvent folder
input_template="$simulation_path/input.txt"
input_file="$solvent_folder/input.txt"
colour_template="$simulation_path/COLOUR.txt"
colour_file="$solvent_folder/COLOUR.txt"

if [[ -f "$input_template" ]]; then
    sed -e "s|FILEPATH|${solvent_folder}/|g" \
        -e "s|SIM_NAME|${solvent}|g" \
        "$input_template" > "$input_file"
    echo "COPYING COLOUR FILE"
    cp $colour_template $colour_file
else
    echo "Template input.txt file not found in the simulation path."
    exit 1
fi

input_path="../addinput_stg_2.txt"

# Run the CrystoGen simulation with the found input file
echo "CrystoGen initiated at $(basename "$solvent_folder")"
cd "$solvent_folder"

"$CG_EXE" < "$input_path" > "${solvent_folder}/stdout.txt" 2> "${solvent_folder}/stderr.txt"
echo "Simulation Complete: $(basename "$solvent_folder")"

echo "DONE!"

