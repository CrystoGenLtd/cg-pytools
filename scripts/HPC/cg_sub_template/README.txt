Template for submitting a single CG job.

Will make a subfolder called "seed" and run a CG job using the input.txt, addinput.txt and net.txt files in the same folder.

Make sure to add a net.txt and STRUCTURE_FILE.txt to the folder for the structure you are modelling.

input.txt and net.txt must keep their names. 
addinput.txt and STRUCTURE_FILE.txt can be changed, but these must be updated in the cg_folder_sub.sh and input.txt, respectively.

Otherwise, just run sbatch cg_folder_sub.sh and it will create a subfolder, populate input values and submit the calculation.

Useful for growing single seeds, or running longer calculations on difficult structures where screening is costly.