Setup for running a batch of CG calculations with different net files.

Reads from a file called "file_list.txt" to find the paths for the net files to use in order.

Requires input.txt, addinput.txt, structure file in directory referenced by input.txt.
Optional - checkpoint and colouring files referenced in input.txt

Submit with sbatch mixed_screen.sh to generate a folder with each net file and run a CG calculation.