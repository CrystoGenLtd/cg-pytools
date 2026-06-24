#!/bin/bash -l
#$ -cwd
#$ -pe smp.pe 4

# Set environment variables, path to CIF, and OCC/CG executables  
CG_KEY="YOUR_PATH/crystogen-1.2.0-linux-x86_64/bin/CrystoGen.key"
CG_EXE="YOUR_PATH/crystogen-1.2.0-linux-x86_64/bin/crystogen"

input_path="../addinput.txt"

# Run the CrystoGen simulation with the found input file
"$CG_EXE" < "$input_path" > "./stdout.txt" 2> "./stderr.txt"

echo "DONE!"

