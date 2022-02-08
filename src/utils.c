#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/resource.h>
#include <errno.h>
#include "utils.h"

void print_error(const char *err) {
    printf("\n\nERROR: %s \n\n", err);
    fflush(NULL);
    exit(1);
}

void print_error_code(const char *err, int status) {
    printf("\n\nERROR: %s -> Error code: %d\n\n", err, status);
    fflush(NULL);
    exit(status);
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
    int pos = i * inst->nnodes + j - ((i + 1) * (i + 2) / 2);
    return pos;
}

int xxpos(int i, int j, instance *inst) {
    int pos = i * (inst->nnodes) + j;
    return pos;
}

int ypos(int i, int j, instance *inst) {
    int xpos = inst->nnodes * inst->nnodes;
    int ypos = xpos + (i * (inst->nnodes) + j);
    return ypos;
}

void print_stats(instance *inst, double time_passed) {
    FILE *stats = fopen("../logfiles/stats.txt", "a");
    if (stats == NULL) {
        print_error_code("Error openining stats files", ENOENT);
    }

    fprintf(stats, "Problem-> %-40s", inst->input_file);
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
    return (int) (r * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
}

double att_distance(int i, int j, instance *inst) {
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
    } else if (strncmp(inst->edge_type, "EUC_2D", 6) == 0) {
        return euc_2d_distance(i, j, inst);
    } else if (strncmp(inst->edge_type, "GEO", 3) == 0) {
        return geo_distance(i, j, inst);
    }
}

double seconds() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) ts.tv_sec + 1.0e-9 * ((double) ts.tv_nsec);
    return ((double) clock() / (double) CLOCKS_PER_SEC);
}

char *logfilename(instance *inst) {
    //Date info for the log of cplex
    time_t rawtime;
    time(&rawtime);
    struct tm *time_struct = localtime(&rawtime);
    char time_string[20];
    strftime(time_string, 20, "%y_%m_%d-%H_%M", time_struct);

    //Saves the info into the dest array
    char *filename = malloc(strlen("../logfiles/cplex_log_") + 2 + strlen(time_string));
    sprintf(filename, "../logfiles/cplex_log%d_%s", inst->model_type, time_string);
    return filename;
}

void free_instance(instance *inst, int free_solution) {
    free(inst->latitude);
    free(inst->longitude);
    free(inst->x_coord);
    free(inst->y_coord);
    free(inst->component);
    free(inst->successors);
    if (free_solution) free(inst->solution);
}

void compute_solution_from_successors(instance *inst, double *x, int *successors) {
    for (int i = 0; i < inst->nvariables; i++) {
        x[i] = 0.0;
    }

    for (int i = 0; i < inst->nnodes; i++) {
        x[xxpos(i, successors[i], inst)] = 1.0;
    }

}

double compute_solution_cost(instance *inst, int *successors) {
    double cost = 0.0;
    for (int i = 0; i < inst->nnodes; i++) {
        cost += distance(i, successors[i], inst);
    }
    return cost;
}

void compute_bigger_cut(instance *inst, int *first_node, int *second_node, double *improvement) {
    *improvement = 0;
    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = 0; j < inst->nnodes; j++) {
            //Skip rules
            if (inst->successors[i] == -1
                || inst->successors[j] == -1
                || i == j
                || inst->successors[i] == j
                || inst->successors[j] == i)
                continue;

            //Difference of the 4 edges taken into considerations
            double delta = distance(i, inst->successors[i], inst) +
                           distance(j, inst->successors[j], inst) -
                           distance(i, j, inst) -
                           distance(inst->successors[i], inst->successors[j], inst);

            //If there is an improvement I save the values of the nodes to be changed
            if (delta > *improvement) {
                *first_node = i;
                *second_node = j;
                *improvement = delta;
            }
        }
    }
}

void perfrom_cut(instance *inst, const int *first_node, const int *second_node, const double *improvement){
    inst->best_value -= *improvement;
    int current = inst->successors[*first_node];
    int previous = inst->successors[*second_node];
    int end_node = previous;
    while (current != end_node) {
        int next = inst->successors[current];
        inst->successors[current] = previous;
        previous = current;
        current = next;
    }
    inst->successors[*first_node] = *second_node;
}