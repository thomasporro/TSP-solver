#define _CRT_SECURE_NO_DEPRECATE
#define PI 3.141592
#define MAX_TIME_LIMIT 3600.0
#define N_ENUMS 18

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cplex.h>
#include <float.h>
#include <errno.h>
#include "read_input.h"
#include "utils.h"

void parse_command_line(int argc, char **argv, instance *inst) {
    //Default values
    inst->timelimit = MAX_TIME_LIMIT;
    inst->model_type = -1;
    inst->performance_profile = 0;
    inst->model_type_counter = 0;
    int flag_perf_prof = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-file") == 0) {
            strcpy(inst->input_file, argv[++i]);
            printf("Input file name: %s\n", inst->input_file);
            flag_perf_prof = 0;
            continue;
        }
        if (strcmp(argv[i], "-model_type") == 0) {
            inst->model_type = atoi(argv[++i]);
            printf("Model type selected: %d\n", inst->model_type);
            flag_perf_prof = 0;
            continue;
        }
        if (strcmp(argv[i], "-time_limit") == 0) {
            inst->timelimit = atoi(argv[++i]);
            printf("Time limit selected: %f\n", inst->timelimit);
            flag_perf_prof = 0;
            continue;
        }
        if (strcmp(argv[i], "-perf_prof") == 0) {
            inst->performance_profile = atoi(argv[++i]);
            inst->model_type_vector = (int *) calloc(N_ENUMS, sizeof(int));
            printf("Performance profile: %d\n", inst->performance_profile);
            flag_perf_prof = 1;
            continue;
        }
        if(flag_perf_prof){
            inst->model_type_vector[inst->model_type_counter++] = atoi(argv[i]);
        }
    }
};

void read_input(instance *inst) {
    //Declaration of local variable
    char *parameter_name;
    char line[180];
    char *token;
    int active_session = 0;

    //Initialization of the values into inst
    inst->nnodes = -1;
    inst->nvariables = -1;
    inst->best_value = DBL_MAX;

    //Opens the files that will be read
    FILE *file = fopen(inst->input_file, "r");
    if (file == NULL) {
        print_error_code("No such input file", ENOENT);
    }

    //Cycles for each line into the input file
    while (fgets(line, sizeof(line), file) != NULL) {
        if (strlen(line) <= 1) { continue; }

        parameter_name = strtok(line, " :");

        //Will save the number of nodes and alloc the memory for the
        //coordinate's value
        if (strncmp(parameter_name, "DIMENSION", 9) == 0) {
            active_session = 0;
            inst->nnodes = atoi(strtok(NULL, " :"));
            inst->x_coord = (double *) calloc(inst->nnodes, sizeof(double));
            inst->y_coord = (double *) calloc(inst->nnodes, sizeof(double));

            //REMOVE
            inst->latitude = (double *) calloc(inst->nnodes, sizeof(double));
            inst->longitude = (double *) calloc(inst->nnodes, sizeof(double));

            //Alloc the memory for successors and component
            inst->successors = (int *) calloc(inst->nnodes, sizeof(int));
            inst->component = (int *) calloc(inst->nnodes, sizeof(int));

            //Print of debug
            printf("Number of nodes: %d\n", inst->nnodes);
        }

        //Saves the edge type into the instance
        if (strncmp(parameter_name, "EDGE_WEIGHT_TYPE", 16) == 0) {
            strcpy(inst->edge_type, strtok(NULL, " :"));

            //Print of debug
            printf("Edge weight type: %s\n", inst->edge_type);
        }

        //Will redirect the parser to the active session so we can parse
        //the nodes coordinates
        if (strncmp(parameter_name, "NODE_COORD_SECTION", 18) == 0) {
            active_session = 1;
            continue;
        }

        //Will end the reading of the file
        if (strncmp(parameter_name, "EOF", 3) == 0) {
            active_session = 0;
            break;
        }

        //Parse the coordinates
        if (active_session == 1) {
            int nodes_number = atoi(parameter_name) - 1;

            //X's coordinates
            token = strtok(NULL, " :,");
            inst->x_coord[nodes_number] = atof(token);

            if (strncmp(inst->edge_type, "GEO", 3) == 0) {
                int deg = round(inst->x_coord[nodes_number]);
                double min = inst->x_coord[nodes_number] - deg;
                inst->latitude[nodes_number] = PI * (deg + 5.0 * min / 3.0) / 180.0;
            }

            //Y's coordinates
            token = strtok(NULL, " :,");
            inst->y_coord[nodes_number] = atof(token);

            if (strncmp(inst->edge_type, "GEO", 3) == 0) {
                int deg = round(inst->y_coord[nodes_number]);
                double min = inst->y_coord[nodes_number] - deg;
                inst->longitude[nodes_number] = PI * (deg + 5.0 * min / 3.0) / 180.0;
            }

            continue;
        }
    }
};
