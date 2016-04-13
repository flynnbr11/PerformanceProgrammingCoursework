This directory contains the files to run a molecular dynamics simulation, 
	developed as course work for Performance Programming in submission for 
	MSc in High Performance Computing, University of Edinburgh.
	
To compile the optimised program: 
	$ make
	
And run via
	$ ./MD

To run on Morar (must compile first)
	$ qsub runMorar.sge	
	
Files contained within: 
* control.c:
		Contains main function. Changed from original as outlined in report. 
* MD.c:
		Contains evolve function. Changed from original as outlined in report. 
* coord.h: 
		Contains initialisations of variables for use throughout simulation. 
		Changed from original as outlined in report. 
* util.c: 
		Contains functions. No explicit purpose since functions are all inlined, 
		but retained for backwards compatability and reference. 
* Test:
		Directory to check if output is correct. 
		To check correctness of most recent run of program:
			$ cd Test
			$ ./diff-output correctOutputs/output.dat200 ../output.dat200
			

