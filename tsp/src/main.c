#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "tsp.h"
#include "read_input.h"
#include "plot.h"
#include "utils.h"

int performance_profile(instance *inst, int *models, int nmodels, double time_limit);


int main(int argc, char **argv) {
	//Variables
	instance inst;

	//If the command line arguments are minor than 2 exits from the program
	if (argc < 2) {
		printf("Usage error");
		exit(1);
	}

	int model_type[] = { MTZ, MTZ_LAZY, MTZ_IND, GG };
	//return performance_profile(&inst, &model_type, 4, 3600.0);


	//Parse the command line and read the input file
	printf("---------------INPUT FILE INFORMATIONS---------------\n");
	parse_command_line(argc, argv, &inst);
	read_input(&inst);
	

	//Calculate the solution of the problem
	printf("\n--------------OPTIMIZATION INFORMATIONS--------------\n");
	TSPopt(&inst);

	
	//Setting the commands to pass to gnuplot to print the graph
	char *commandsForGnuplot[3];
	commandsForGnuplot[0] = "set title \"GRAPH\"";
	if (inst.model_type == STANDARD || inst.model_type == BENDERS || inst.model_type == BRANCH_AND_CUT) {
		commandsForGnuplot[1] = "plot \"data.dat\" with linespoints linestyle 1 lc rgb \"red\"";
	}
	else {
		commandsForGnuplot[1] = "plot \"data.dat\" using 1:2 with points ls 5 lc rgb \"red\", \\\n"
			"\"arrows.dat\" using 1:2:3:4 with vectors filled head lc rgb \"black\"";
	}

	//Plot the solution with the passed commands
	printf("\n----------------------PLOTTING------------------------\n");
	plot(commandsForGnuplot, 2, &inst);
	
	
	return 0;
}


int performance_profile(instance *inst, int *models, int nmodels, double time_limit) {
	printf("------------START PERFORMANCE PROFILE MODE-----------\n");
	double start_time = seconds();
	inst->timelimit = time_limit;
	printf("TIME LIMIT SETTED TO: %6.2fs\n", inst->timelimit);

	FILE *csv = fopen("performance_profile.csv", "w");
	fprintf(csv, "%d, ", nmodels);

	//Print the first line of the csv file
	for (int i = 0; i < nmodels; i++) {
		fprintf(csv, "%d", models[i]);
		if (i != nmodels - 1)
			fprintf(csv, ", ");
		else
			fprintf(csv, "\n");
	}

	FILE *files = fopen("files.txt", "r");
	char line[180];
	char *file_name;
	while (fgets(line, sizeof(line), files) != NULL) {
		if (strlen(line) <= 1) { continue; }

		//Retrieve the file name removing the newline
		file_name = strtok(line, "\n");

		printf("\nTESTING FILE: %s\n", file_name);
		fprintf(csv, "%s, ", file_name);

		//Read and parse the file to test
		strcpy(inst->input_file, file_name);
		read_input(inst);

		//Iterate over the models passed
		for (int i = 0; i < nmodels; i++) {
			printf("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n");
			inst->model_type = models[i];
			double start_time = seconds();
			TSPopt(inst);
			double end_time = seconds();

			//Print
			fprintf(csv, "%f", end_time - start_time);
			if (i != nmodels - 1)
				fprintf(csv, ", ");
			else
				fprintf(csv, "\n");
			printf("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n\n");
		}
	}

	double end_time = seconds();
	printf("TIME ELAPSED: %f s\n\n", end_time - start_time);

	return 0;
}