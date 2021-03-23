#include <stdlib.h>
#include <stdio.h>
#include "tsp.h"
#include "read_input.h"
#include "plot.h"

int main(int argc, char **argv) {
	//Variables
	instance inst;

	//If the command line arguments are minor than 2 exits from the program
	if (argc < 2) {
		printf("Usage error");
		exit(1);
	}

	//Parse the command line and read the input file
	printf("---------------INPUT FILE INFORMATIONS---------------\n");
	parse_command_line(argc, argv, &inst);
	read_input(&inst);

	int model_type[] = {0, 10, 11, 12, 20};

	//Calculate the solution of the problem
	printf("\n--------------OPTIMIZATION INFORMATIONS--------------\n");
	TSPopt(&inst);

	/*
	for (int i = 0; i < 5; i++) {
		inst.model_type = model_type[i];
		TSPopt(&inst);
	}*/

	//Setting the commands to pass to gnuplot to print the graph
	char *commandsForGnuplot[3];
	commandsForGnuplot[0] = "set title \"GRAPH\"";
	if (inst.model_type == 0) {
		commandsForGnuplot[1] = "plot \"data.dat\" with linespoints linestyle 1 lc rgb \"red\"";
	}
	else {
		commandsForGnuplot[1] = "plot \"data.dat\" using 1:2 with points ls 5 lc rgb \"red\", \\\n"
			"\"arcs.dat\" using 1:2:3:4 with vectors filled head lc rgb \"black\"";
	}

	//Plot the solution with the passed commands
	printf("\n----------------------PLOTTING------------------------\n");
	plot(commandsForGnuplot, 2, &inst);

	
	return 0;
}