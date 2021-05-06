#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
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
        print_error("Usage error");
    }


    //Parse the command line and read the input file
    printf("---------------INPUT FILE INFORMATIONS---------------\n");
    parse_command_line(argc, argv, &inst);

    if(inst.performance_profile){
        int model_type[] = {BENDERS, BRANCH_AND_CUT, BRANCH_AND_CUT_RLX};
        return performance_profile(&inst, (int *) &model_type, 3, 3600.0);
    }

    read_input(&inst);


    //Calculate the solution of the problem
    printf("\n--------------OPTIMIZATION INFORMATIONS--------------\n");
    TSPopt(&inst);


    //Setting the commands to pass to gnuplot to print the graph
    char *commandsForGnuplot[3];
    commandsForGnuplot[0] = "set title \"GRAPH\"";
    if (inst.model_type == STANDARD || inst.model_type == BENDERS || inst.model_type == BRANCH_AND_CUT ||
        inst.model_type == DEFAULT || inst.model_type == HARD_FIX_BAC || inst.model_type == SOFT_FIX
        || inst.model_type == BRANCH_AND_CUT_RLX) {
        commandsForGnuplot[1] = "plot \"../testfiles/data.dat\" with linespoints linestyle 1 lc rgb \"red\"";
    } else {
        commandsForGnuplot[1] = "plot \"../testfiles/data.dat\" using 1:2 with points ls 5 lc rgb \"red\", \\\n"
                                "\"../testfiles/arrows.dat\" using 1:2:3:4 with vectors filled head lc rgb \"black\"";
    }

    //Plot the solution with the passed commands
    printf("\n----------------------PLOTTING------------------------\n");
    plot(commandsForGnuplot, 2, &inst);

    free_instance(&inst);
    return 0;
}


int performance_profile(instance *inst, int *models, int nmodels, double time_limit) {
    printf("------------START PERFORMANCE PROFILE MODE-----------\n");

    //Setting parameters inside of inst
    inst->performance_profile = 1;
    inst->timelimit = time_limit;
    printf("TIME LIMIT SETTED TO: %6.2fs\n", inst->timelimit);
    double start_time = seconds();

    FILE *csv = fopen("../logfiles/performance_profile.csv", "w");
    if (csv == NULL) {
        print_error_code("Failed to open ../logfiles/performance_profile.csv", ENOENT);
    }
    fprintf(csv, "%d, ", nmodels);

    //Print the first line of the csv file
    for (int i = 0; i < nmodels; i++) {
        fprintf(csv, "%d", models[i]);
        if (i != nmodels - 1)
            fprintf(csv, ", ");
        else
            fprintf(csv, "\n");
    }

    FILE *files = fopen("../testfiles/files.txt", "r");
    if (files == NULL) {
        print_error_code("Failed to open ../testfiles/files.txt", ENOENT);
    }

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
            double start_time_for = seconds();
            TSPopt(inst);
            double end_time_for = seconds();

            //Print
            fprintf(csv, "%f", end_time_for - start_time_for);
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