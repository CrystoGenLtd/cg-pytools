This is a script to run a single CG job in a folder.

To run a CG job in a folder, you need:
input.txt - this is the main input for CG which has a strict format and must be called input.txt
net.txt - this is the energy data for interactions in the structure which also must be called net.txt

addinput.txt - this is a variable file which contains data in response to some of the prompts for the settings
in input.txt (e.g. if you choose delta mu mode 3 in input.txt, then you need to tell CG how long to stay at
high driving force, and how long to take when equilibrating - these come from addinput.txt). 
The script here relies on this being called addinput.txt and in the same folder, but this can be changed.

STRUCTURE_FILE.txt - this is called from input.txt, so you can change this path to whatever works best for you.

Completed GUI simulations write out input and addinput files for you to make things easier.

Optional:
checkpoint files, if you're using a checkpoint simulation, will need to be provided and the path defined in input.txt
colouring files, the same is required here if using colouring files or looking for size output.
