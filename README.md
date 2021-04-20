# TSP-solver

This project is the result of the course of Linear Programming 2 @unipd. The main goal of this program is to solve the Traveling Salesman Problem using the
powerful CPLEX by IBM.


This code comes with a makefile in order to compile it by yourself.

## Usage
`-f <namefile>` This pass to the program the file you want to use. It's supported only the *.tsp models with nodelist. 

`-time_limit <timelimit>` Sets the timelimit to the program.

`-model_type <modeltype>` This pass to the program the model you want to use during the optimization. Here's a list of the supported ones:
- STANDARD = 0
- BENDERS = 1
- BRANCH AND CUT = 2
- MTZ = 10
- MTZ WITH LAZY CONSTRAINTS = 11
- MTZ WITH INDICATOR CONSTRAINTS = 12
- GG = 20
