#!/bin/bash -l
#SBATCH -p QUEUE_NAME
#SBATCH -c NUM_CORES

# Use a default input deck unless specified on the command-line
if [ -n "$1" ]; then
  CIF=$1
else
  # Need to be in the current dir
  CIF=MY_CIF.cif
fi

# Activate python environment if using cg_net.py to group interactions
source PYTHON_ENV_PATH

# Set environment variables, path to CIF, and OCC/CG executables  
CG_KEY="KEY_PATH"
CG_EXE="CG_EXECUTABLE_PATH"

# Extract Solvent based on JobArrayID
SOLVENT="water"

SIM_DIR=$(pwd)
cif_base=$(basename "$CIF" .cif)

# echo "SOLVENT: ${SOLVENT}"
# echo "SIM_DIR: ${SIM_DIR}"

# Path to the directory where the simulation should run
solvent="${SOLVENT// /+}"
simulation_path="$SIM_DIR"

# Convert to absolute path
simulation_path=$(realpath "$SIM_DIR")
input_path="$simulation_path/addinput.txt"

# Check if input file was found
if [[ ! -f "$input_path" ]]; then
    echo "Additional input file not found in parent directories."
    exit 1
fi

# Create solvent-specific folder
solvent_folder="$simulation_path/seed"
mkdir -p "$solvent_folder"

# Copy the solvent_net.txt file to solvent folder as net.txt
# Option1: If no grouping is required:
# cp "$simulation_path/${solvent}_net.txt" "$solvent_folder/net.txt"
# Option2: If interactions need grouping:
python cg_net.py -i "$simulation_path/net.txt" -o "$solvent_folder/net.txt" --group

ln -s "$CG_KEY" "$solvent_folder/CrystoGen.key"

# Modify and copy input.txt file to solvent folder
input_template="$simulation_path/input.txt"
input_file="$solvent_folder/input.txt"

if [[ -f "$input_template" ]]; then
    sed -e "s|FILEPATH|${solvent_folder}/|g" \
        -e "s|SIM_NAME|${solvent}|g" \
        "$input_template" > "$input_file"
else
    echo "Template input.txt file not found in the simulation path."
    exit 1
fi

# Run the CrystoGen simulation with the found input file
echo "CrystoGen initiated at $(basename "$solvent_folder")"
cd "$solvent_folder"
"$CG_EXE" < "$input_path" > "${solvent_folder}/stdout.txt" 2> "${solvent_folder}/stderr.txt"
echo "Simulation Complete: $(basename "$solvent_folder")"

echo "DONE!"

echo "CREATING TEMPLATE COLOUR FILE"

colour=$(find . -type f -name "*all.txt" -not -empty -print -quit)
lines=$(wc -l < $colour)

if [ -f ../COLOUR.txt ]; then
    rm ../COLOUR.txt
fi

echo " ${lines}" >> ../COLOUR.txt
cat $colour >> ../COLOUR.txt
