NOTE - NEW SCRIPTS RELY ON NEW INSTALLATION OF OCC
OCC can now be installed as a python package (occ-py) - removes need to provide paths for occ files.

Two stage screening process to ensure a facets-all file is produced. No longer needs an internal licence key - seed engineering is available for all users.

All three scripts have changed slightly, so make sure to copy all three and change your paths / python environment.

job.sh - runs the initial OCC calculation for the solid phase 
- copies and renames files to STRUCTURE_FILE.txt
- copies generated net file to a template file called net.txt

screen.sh - runs a single calculation with the generated net file to try and force a facet-all file to be written.
- short calculation so it runs very quickly, and limited to facets from -2,-2,-2 to 2,2,2.
- runs in the "seed" folder by default.
- copies the template COLOUR.txt file to the parent folder to be used by all following calculations.
- takes input.txt and addinput.txt as inputs for the template CG calculation.
- new addition uses the crystal cutting feature to restrict the crystal to grow within a sphere of 10 unit cells.

screen_stg_2.sh - runs the complete screen (179 solvents) with the COLOUR.txt file as an input.
- takes input_stg_2.txt and addinput_stg_2.txt as inputs for the growth calculations.
- outputs a size measurement to the size.csv file every 1% of the simulation (e.g. every 30k steps if 3 million iterations are set).


submit_all.sh will submit everything in order as previous.

run with "sbatch submit_all.sh CIF_FILE.cif"

Troubleshooting:
Box size issues, CG should catch if the box edge is reached, but if this is becoming a problem for NON-NEEDLELIKE crystals, then the box can be adjusted manually. Set the automatic box size in input_stg_2.txt to "no" then add unit cell size for a, b, and c after the third line in addinput_stg_2.txt (e.g.
yes
no
801
801
801
100000 etc.)
