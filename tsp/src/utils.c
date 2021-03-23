#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include "utils.h"

void print_error(const char *err) {
	printf("\n\nERROR: %s \n\n", err);
	fflush(NULL);
	exit(1);
}

double distance(int i, int j, instance *inst) {
	//Differences between the x's and y's coordinates
	double x_dist = inst->x_coord[i] - inst->x_coord[j];
	double y_dist = inst->y_coord[i] - inst->y_coord[j];

	//Calculates the Euclidian Distance and returns it
	int value = sqrt(x_dist * x_dist + y_dist * y_dist) + 0.499999999;
	return value + 0.0;
}

int xpos(int i, int j, instance *inst) {
	if (i == j) print_error("i == j in xpos");
	if (i > j) return xpos(j, i, inst);
	int pos = i * inst->nnodes + j - ((i + 1)*(i + 2) / 2);
	return pos;
}

int xxpos(int i, int j, instance *inst) {
	int pos = i * (inst->nnodes) + j;
	return pos;
}

int ypos(int i, int j, instance *inst) {
	int xpos = inst->nnodes  * inst->nnodes;
	int ypos = xpos + (i * (inst->nnodes) + j);
	return ypos;
}

void print_stats(instance *inst, double time_passed) {
	FILE *stats = fopen("stats.txt", "a");

	fprintf(stats, "Problem -> %s\n", inst->input_file);
	fprintf(stats, "Model -> %d\n", inst->model_type);
	fprintf(stats, "Time passed -> %f\n\n", time_passed);

	fclose(stats);
}