These scripts run a solvent screen with OCC and CG
Requires a CIF to start

Run job.sh with sbatch job.sh CIF_NAME.cif

This will run a single OCC calculation in water to compute the lattice energies (solid phase)
This should be run to completion first, as OCC will then be able to reuse the lattice energy for all solvent calculations without repeating.

Once finished, a folder containing dimers will be generated, you can visualise these as XYZ files if required.

Then run screen.sh with sbatch screen.sh CIF_NAME.cif

This will run through the solvent.txt folder, compute the net files for CG with OCC in each solvent.
Then it will run a single CG calculation with said energies to grow a single crystal.

note: requires a virtual python or conda environment to run. See pdf manual for details.