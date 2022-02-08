# TSP-solver

This project is the result of the course of Linear Programming 2 @unipd. The main goal of this program is to solve the Travelling Salesman Problem.

We implemented the code in order to solve the problem in different ways. We start the course using CPLEX and various model implemented in it. Later on we developed heuristics and meta-heuristic on order to solve really large problem without search for optimality


## Installation

This code comes with a makefile in order to compile it by yourself.

## Usage

In this section you will find the commands you can pass to the program in order to execute the code of a TSP problem. Here I will put some usage example from linux.

The program requires 2 arguments as input, otherwise it will raise an error. Here's the list of commands you can pass to it:

```bash
./tsp [-file <namefile>] [-time_limit <timelimit>] [-model_type <modeltype>] [-perf_prof <flag>]
```


`-file <namefile>` 

This pass to the program the file you want to use. It's supported only the *.tsp models with explicit nodelist. 

`-time_limit <timelimit>` 

Sets the timelimit to the program. You must pass the timelimit in seconds.

`-model_type <modeltype>` 

This pass to the program the model you want to use during the optimization. Here's a list of the supported ones:
- standard = 0;
- bender's method = 1;
- branch and cut = 2;
- branch and cut + relaxation = 3;
- hard fixing = 4;
- soft fixing = 5;
- mtz = 10;
- mtz with lazy constraints = 11;
- mtz with indicator constraints = 12;
- gg = 20;
- greedy = 30;
- greedy with 2-opt refining = 31;
- extra mileage = 32;
- vns = 40.

## Report (in continuous update)
You can download the report with this [link](https://github.com/thomasporro/TSP-solver/blob/greedy/latex_report/ro2.pdf).
