#define _CRT_SECURE_NO_DEPRECATE
#define EPS 1e-5

#include "plot.h"
#include "utils.h"

void plot(char **commands, int n_commands, instance *inst) {
	//Open the gnuplot enviroment and write the data into a data
	FILE *gnuplotPipe = _popen("C:/UNIPD/ro2/gnuplot/bin/gnuplot.exe -persistent", "w");
	FILE *temp = fopen("data.dat", "w");
	FILE *arcs = fopen("arrows.dat", "w");

	switch (inst->model_type)
	{
	case 0:
		print_st(temp, inst);
		break;
	case 1:
		print_st(temp, inst);
		break;
	case 10:
		print_MTZ(temp, arcs, inst);
		break;
	case 11:
		print_MTZ(temp, arcs, inst);
		break;
	case 12:
		print_MTZ(temp, arcs, inst);
		break;
	case 20:
		print_MTZ(temp, arcs, inst);
		break;
	default:
		print_st(temp, inst);
		break;
	}

	//Closing the file
	fclose(temp);
	fclose(arcs);

	//Executing the commands passed into the function
	for (int i = 0; i < n_commands; i++) {
		fprintf(gnuplotPipe, "%s \n", commands[i]);
	}

	//Closing the pipeline of gnuplot
	_pclose(gnuplotPipe);
}

void print_st(FILE *temp, instance *inst) {
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = i + 1; j < inst->nnodes; j++) {
			int k = xpos(i, j, inst);

			if (fabs(inst->solution[k]) > EPS && fabs(inst->solution[k] - 1.0) > EPS)
				print_error("wrong solution in print_st()");

			//let's modify the value of the edges in each node
			if (inst->solution[k] > 0.5) {
				fprintf(temp, "%lf %lf \n", inst->x_coord[i], inst->y_coord[i]);
				fprintf(temp, "%lf %lf \n", inst->x_coord[j], inst->y_coord[j]);
				fprintf(temp, "\n");
			}
		}
	}
}

void print_MTZ(FILE *temp, FILE *arcs, instance *inst) {
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = 0; j < inst->nnodes; j++) {
			if (i == j) continue;
			int k = xxpos(i, j, inst);

			if (fabs(inst->solution[k]) > EPS && fabs(inst->solution[k] - 1.0) > EPS)
				print_error("wrong solution in print_st()");

			//let's modify the value of the edges in each node
			if (inst->solution[k] > 0.5) {
				//Points elements
				fprintf(temp, "%lf %lf \n", inst->x_coord[i], inst->y_coord[i]);
				fprintf(temp, "%lf %lf \n", inst->x_coord[j], inst->y_coord[j]);
				fprintf(temp, "\n");

				//Vectors element (x, y) and (delta_x, delta_y)
				fprintf(arcs, "%lf %lf ", inst->x_coord[i], inst->y_coord[i]);
				fprintf(arcs, "%lf %lf \n", inst->x_coord[j] - inst->x_coord[i], inst->y_coord[j] - inst->y_coord[i]);
			}
		}
	}
}