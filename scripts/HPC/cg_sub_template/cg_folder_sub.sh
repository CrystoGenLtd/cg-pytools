#!/bin/bash -l
#SBATCH -p QUEUE_NAME
#SBATCH -c NUM_CORES

# Set environment variables, path to CIF, and OCC/CG executables  
CG_KEY="CG_LICENCE_KEY_PATH"
CG_EXE="CG_EXECUTABLE_PATH"

SIM_DIR=$(pwd)
simulation_path="$SIM_DIR"

# Convert to absolute path
simulation_path=$(realpath "$SIM_DIR")
input_path="$simulation_path/addinput.txt"

# Create solvent-specific folder
solvent="seed"
solvent_folder="$simulation_path/seed"
mkdir -p "$solvent_folder"

# copy net file
cp net.txt $solvent_folder

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

