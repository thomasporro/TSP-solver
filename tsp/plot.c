#define _CRT_SECURE_NO_DEPRECATE

#include "plot.h"

void plot(char **commands, int n_commands, instance *inst) {
	//Open the gnuplot enviroment and write the data into a data
	FILE *gnuplotPipe = _popen("C:/UNIPD/ro2/gnuplot/bin/gnuplot.exe -persistent", "w");
	FILE * temp = fopen("data.dat", "w");

	for (int i = 0; i < inst->nnodes; i++) {
		fprintf(temp, "%lf %lf \n", inst->x_coord[i], inst->y_coord[i]);
		fprintf(temp, "e\n");
	}

	//Closing the file
	fclose(temp);

	//Executing the commands passed into the function
	for (int i = 0; i < n_commands; i++) {
		fprintf(gnuplotPipe, "%s \n", commands[i]);
	}

	//Closing the pipeline of gnuplot
	_pclose(gnuplotPipe);
}