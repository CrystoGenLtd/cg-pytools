NOTE - NEW SCRIPTS RELY ON NEW INSTALLATION OF OCC
OCC can now be installed as a python package (occ-py) - removes need to provide paths for occ files.

New growth rate scripts - a bit more flexible
Runs a series of CrystoGen calculations starting from the same seed but held at a constant supersaturation.
This simulation will also create a summary.csv file required to plot the the results by supersaturation. 

Submit with ./run_growth.sh starting_delta_mu(kcal/mol) increment_value(kcal/mol) number_of_simulations_to_queue

e.g. ./run_growth.sh -8 0.5 40

Will queue 40 simulations running from -8 kcal/mol to -8 + (0.5*39) = 11.5 kcal/mol

You will need to have growth_rate_screen.sh, input.txt and addinput.txt in the same folder.
You will also need COLOUR.txt (file containing planes to measure growth rates) and CHECKPOINT.txt (a seed crystal to start from).
You will need STRUCTURE_FILE.txt (structure file for the system).

To perform further analysis with the screening scripts, it is also useful to copy in a xx_cg_results.json file from a solvent screen.


Possible issues:
Make sure the run_growth.sh is an executable script (chmod u+x run_growth.sh)
Make sure that the growth_rate_screen.sh script parameters match your HPC system.
If the simulation box edge is continually being reached - try adjusting the box size in addinput.txt (e.g. 801 instead of 401), this will require more memory, but will allow larger crystals to grow.