#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "utils.h"

void print_error(const char *err) {
	printf("\n\nERROR: %s \n\n", err);
	fflush(NULL);
	exit(1);
}

double euc_2d_distance(int i, int j, instance *inst) {
	//Differences between the x's and y's coordinates
	double x_dist = inst->x_coord[i] - inst->x_coord[j];
	double y_dist = inst->y_coord[i] - inst->y_coord[j];

	//Calculates the Euclidian Distance and returns it
	int value = round(sqrt(pow(x_dist, 2) + pow(y_dist, 2)));
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

	fprintf(stats, "Problem-> %-30s", inst->input_file);
	fprintf(stats, "Model-> %-10d", inst->model_type);
	fprintf(stats, "Time passed-> %-10f\n", time_passed);

	fclose(stats);
}

double geo_distance(int i, int j, instance *inst) {
	double r = 6378.388;
	double pi = 3.141592;

	
	double q1 = cos(inst->longitude[i] - inst->longitude[j]);
	double q2 = cos(inst->latitude[i] - inst->latitude[j]);
	double q3 = cos(inst->latitude[i] + inst->latitude[j]);
	return (int)(r * acos(0.5*((1.0 + q1)*q2 - (1.0 - q1)*q3)) + 1.0);

	/*
	double delta_lat = (inst->x_coord[j] - inst->x_coord[i]) * PI / 180.0;
	double delta_lon = (inst->y_coord[j] - inst->y_coord[i]) * PI / 180.0;

	double start_lat = inst->x_coord[i] * PI / 180.0;
	double end_lat = inst->x_coord[j] * PI / 180.0;

	double formula = pow(sin(delta_lat / 2), 2) + pow(sin(delta_lon / 2), 2)*cos(start_lat)*cos(end_lat);
	double c = 2 * asin(formula);
	return c * R;*/
}

double att_distance(int i, int j, instance *inst){
	double delta_x = inst->x_coord[i] - inst->x_coord[j];
	double delta_y = inst->y_coord[i] - inst->y_coord[j];
	double rij = sqrt((pow(delta_x, 2) + pow(delta_y, 2)) / 10.0);
	double tij = round(rij);
	if (tij < rij) return tij + 1;
	else return tij;

}

double distance(int i, int j, instance *inst) {
	if (strncmp(inst->edge_type, "ATT", 3) == 0) {
		return att_distance(i, j, inst);
	}
	else if (strncmp(inst->edge_type, "EUC_2D", 6) == 0) {
		return euc_2d_distance(i, j, inst);
	}
	else if (strncmp(inst->edge_type, "GEO", 3) == 0) {
		return geo_distance(i, j, inst);
	}
}

double seconds(){
	return ((double)clock() / (double)CLOCKS_PER_SEC);
}