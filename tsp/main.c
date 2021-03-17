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

	//Commands for gnuplot
	char *commandsForGnuplot[] = { "set title \"POINTS\"", "plot \"data.dat\"" };

	//Parse the command line and read the input file
	printf("---------------INPUT FILE INFORMATIONS---------------\n");
	parse_command_line(argc, argv, &inst);
	read_input(&inst);

	//Calculate the solution of the problem
	printf("\n--------------OPTIMIZATION INFORMATIONS--------------\n");
	TSPopt(&inst);

	//Plot the solution with the passed commands
	printf("\n----------------------PLOTTING------------------------\n");
	plot(commandsForGnuplot, 2, &inst);

	
	return 0;
}