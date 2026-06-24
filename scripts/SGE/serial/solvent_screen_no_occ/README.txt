These scripts run a solvent screen only using CG
Requires a CIF and net files to have been calculated by OCC already

Run sbatch_no_occ screen with sbatch screen_no_occ.sh CIF_NAME.cif

This will run through the solvent.txt folder, use existing net files from OCC, and make a subfolder for each solvent.
Then it will run a single CG calculation with said energies to grow a single crystal.

note: requires a virtual python or conda environment to run. See pdf manual for details.

Can submit job.sh then screen_no_occ.sh in sequence by submitting submit.sh instead