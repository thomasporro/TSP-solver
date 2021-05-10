#define _CRT_SECURE_NO_DEPRECATE
#define EPS 1e-5

#include <errno.h>
#include "plot.h"
#include "utils.h"

void plot(char **commands, int n_commands, instance *inst) {
    //Open the gnuplot enviroment and write the data into a data
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    FILE *data = fopen("../testfiles/data.dat", "w");
    FILE *arcs = fopen("../testfiles/arrows.dat", "w");
    if (data == NULL || arcs == NULL) {
        print_error_code("Error while opening files in plot()", ENOENT);
    }

    switch (inst->model_type) {
        case STANDARD:
            print_st(data, inst);
            break;
        case BENDERS:
            print_st(data, inst);
            break;
        case BRANCH_AND_CUT:
            print_st(data, inst);
            break;
        case MTZ:
            print_MTZ(data, arcs, inst);
            break;
        case MTZ_LAZY:
            print_MTZ(data, arcs, inst);
            break;
        case MTZ_IND:
            print_MTZ(data, arcs, inst);
            break;
        case GREEDY:
            print_heur(data, inst);
            break;
        case GREEDY_REF:
            print_heur(data, inst);
            break;
        case XTRA_MILEAGE:
            print_heur(data, inst);
            break;
        default:
            print_st(data, inst);
            break;
    }

    //Closing the file
    fclose(data);
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

void print_heur(FILE *temp, instance *inst){
    int start_node;
    for(int i= 0; i<inst->nnodes; i++){
        if(inst->successors[i]!=-1){
            start_node = i;
            break;
        }
    }
    //int start_node = 0;
    int current_node = inst->successors[start_node];

    fprintf(temp, "%lf %lf \n", inst->x_coord[start_node], inst->y_coord[start_node]);
    fprintf(temp, "%lf %lf \n", inst->x_coord[current_node], inst->y_coord[current_node]);
    fprintf(temp, "\n");

    while(start_node != current_node){
        int next_node = inst->successors[current_node];
        fprintf(temp, "%lf %lf \n", inst->x_coord[current_node], inst->y_coord[current_node]);
        fprintf(temp, "%lf %lf \n", inst->x_coord[next_node], inst->y_coord[next_node]);
        fprintf(temp, "\n");
        current_node = next_node;
    }
}