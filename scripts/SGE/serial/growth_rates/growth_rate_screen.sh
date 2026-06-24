#!/bin/bash -l
#$ -cwd

var=$1
inc=$2

# Set environment variables, path to CIF, and OCC/CG executables  
CG_KEY="YOUR_PATH/crystogen-1.2.0-linux-x86_64/bin/CrystoGen.key"
CG_EXE="YOUR_PATH/crystogen-1.2.0-linux-x86_64/bin/crystogen"

# Extract Solvent based on JobArrayID
TASK_ID=${SGE_TASK_ID}

SIM_DIR=$(pwd)

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
job_folder="$simulation_path/$TASK_ID"
mkdir -p "$job_folder"

cp net.txt $job_folder/net.txt

ln -s "$CG_KEY" "$job_folder/CrystoGen.key"

# Modify and copy input.txt file to solvent folder
input_template="$simulation_path/input.txt"
input_file="$job_folder/input.txt"

if [[ -f "$input_template" ]]; then
    var=$(echo "$var + $inc * ($TASK_ID - 1)" | bc)
    sed -e "s|FILEPATH|${job_folder}/|g" \
        -e "s|SIM_NAME|${TASK_ID}|g" \
	-e "s|supersat_vary|${var}|g" \
        "$input_template" > "$input_file";
else
    echo "Template input.txt file not found in the simulation path."
    exit 1
fi

# Run the CrystoGen simulation with the found input file
echo "CrystoGen initiated at $(basename "$job_folder")"
cd "$job_folder"
echo "SUPERSATURATION!"
echo $var
"$CG_EXE" < "$input_path" > "${job_folder}/stdout.txt" 2> "${job_folder}/stderr.txt"
echo "Simulation Complete: $(basename "$job_folder")"

echo "DONE!"

