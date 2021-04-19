#define _CRT_SECURE_NO_DEPRECATE
#define EPS 1e-5

#include "../include/plot.h"
#include "../include/utils.h"

void plot(char **commands, int n_commands, instance *inst) {
	//Open the gnuplot enviroment and write the data into a data
	FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
	FILE *temp = fopen("data.dat", "w");
	FILE *arcs = fopen("arrows.dat", "w");

	switch (inst->model_type)
	{
	case STANDARD:
		print_st(temp, inst);
		break;
	case BENDERS:
		print_st(temp, inst);
		break;
	case BRANCH_AND_CUT:
		print_st(temp, inst);
		break;
	case MTZ:
		print_MTZ(temp, arcs, inst);
		break;
	case MTZ_LAZY:
		print_MTZ(temp, arcs, inst);
		break;
	case MTZ_IND:
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
	pclose(gnuplotPipe);
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