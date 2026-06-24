Run with:
sbatch growth_rate_screen.sh x y

where x = starting supersaturation (kcal/mol)
and y = increment amount (kcal/mol)

Requires the following:

input.txt template
addinput.txt template


STRUCTURE_FILE.txt - your structure file from a screen
net.txt - a net.txt file for your chosen solvent from a screen
COLOUR.txt - file with your facets (e.g. a facets_found.txt file from a previous calc).
CHECKPOINT.txt - a checkpoint file output from a previous calculation to start from.

addinput.txt requires editing - the 401 401 401 box is a rough size estimate for a block crystal.
You can check the checkpoint file for the dimensions that the simulation finished at and scale accordingly.

This line:
#SBATCH --array=1-20

Sets how many simulations you want to run. So this will run from (0 x increment)+starting to (19 x increment)+starting delta mu.

