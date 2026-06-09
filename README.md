# AB-OCEL

This git contains the source code and the evaluation for the object-centric event log sampling algorithm RP-OCEL.

For the evaluation we modified OCPA for multiprocessing and supply the modified code in this git. The original version of OCPA can be found under https://github.com/ocpm/ocpa.

# Running the code

To run the python code for RP-OCEL sampling your need to install python 3.11. A requirements.txt is in the top folder and can be used to install all need packages. 

For sampling:
 - Download the data set from https://zenodo.org/records/13879980 and place the “02_p2p.xml” in the data folder.
 - Start the main.py.
 - When the sample finished successfully the result is in the “output” folder.
 - 
For evaluation:
 - Download the data set from https://zenodo.org/records/13879980 and place the “02_p2p.xml” in the data folder.
 - 	In the main.py uncomment line 16 and 19 to enable the evaluation.
 - Start the main.py
 - The results will be saved in the output folder.
