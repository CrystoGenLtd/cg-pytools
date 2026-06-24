#!/usr/bin/env bash

# Use a default input deck unless specified on the command-line
if [ -n "$1" ]; then
  CIF=$1
else
  # Need to be in the current dir
    CIF=MY_CIF.cif
fi

jobid=$(sbatch --parsable job.sh $CIF)
echo "Submitted initial OCC calculation at ${jobid}"
jobid_2=$(sbatch --dependency=afterok:${jobid} --parsable screen.sh $CIF)
echo "Submitted solvent screening job array!"
sbatch --dependency=afterok:${jobid_2} --parsable screen_stg_2.sh $CIF
echo "Submitted solvent size job array!"
