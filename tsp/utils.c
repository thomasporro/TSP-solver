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
	if (i == j) print_error(" i == j in xpos");
	if (i > j) return xpos(j, i, inst);
	int pos = i * inst->nnodes + j - ((i + 1)*(i + 2) / 2);
	return pos;
}