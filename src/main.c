#define PERFORMANCE_PROFILE 1237030

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "tsp.h"
#include "read_input.h"
#include "plot.h"
#include "utils.h"

/**
 * Function that execute the performance profile with the given information. It saves the results on the
 * performance_profile.csv
 * @param inst A pointer to an instance of the TSP problem (it will be modified)
 * @param models An array containing the
 * @param nmodels The length of the array @param models
 * @param time_limit The time limit of the single run of CPLEX in seconds
 */
void performance_profile(instance *inst, int *models, int nmodels, double time_limit);

/*!
 * Switch to solve the problem with CPLEX or heuristics
 * @param inst is a pointer to the instance where is stored the problem
 * @return 0 if the solution is found. Other values otherwise
 */
int solve(instance *inst);

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

    if (inst.performance_profile) {
        performance_profile(&inst, inst.model_type_vector, inst.model_type_counter, inst.timelimit);
        return 0;
    }

    read_input(&inst);

    //Calculate the solution of the problem
    printf("\n--------------OPTIMIZATION INFORMATIONS--------------\n");
    solve(&inst);

    //Setting the commands to pass to gnuplot to print the graph
    int flag_gnuplot = inst.model_type == STANDARD
                       || inst.model_type == BENDERS
                       || inst.model_type == BRANCH_AND_CUT
                       || inst.model_type == DEFAULT
                       || inst.model_type == HARD_FIX_BAC
                       || inst.model_type == SOFT_FIX
                       || inst.model_type == BRANCH_AND_CUT_RLX
                       || inst.model_type == GREEDY
                       || inst.model_type == GREEDY_REF
                       || inst.model_type == XTRA_MILEAGE
                       || inst.model_type == XTRA_MILEAGE_REF
                       || inst.model_type == VNS
                       || inst.model_type == TABU_SEARCH
                       || inst.model_type == GENETIC;

    char *commandsForGnuplot[3];
    commandsForGnuplot[0] = ""; // Change to change the title of the plot;
    if (flag_gnuplot) {
        commandsForGnuplot[1] = "plot \"../testfiles/data.dat\" with linespoints linestyle 1 lc rgb \"red\"";
    } else {
        commandsForGnuplot[1] = "plot \"../testfiles/data.dat\" using 1:2 with points ls 5 lc rgb \"red\", \\\n"
                                "\"../testfiles/arrows.dat\" using 1:2:3:4 with vectors filled head lc rgb \"black\"";
    }

    //Plot the solution with the passed commands
    printf("\n----------------------PLOTTING------------------------\n");
    plot(commandsForGnuplot, 2, &inst);

    free_instance(&inst, inst.model_type); // !flag_free_solution
    return 0;
}

int solve(instance *inst) {
    double start_time = seconds();
    switch (inst->model_type) {
        case GREEDY:
            printf("Model type chosen: undirected complete graph solved with greedy method\n");
            inst->start_time = seconds();
            greedy(inst);
            printf("Solution cost: %f\n", inst->best_value);
            break;
        case GREEDY_REF:
            printf("Model type chosen: undirected complete graph solved with greedy method + 2-opt refining\n");
            inst->start_time = seconds();
            greedy(inst);
            printf("Greedy cost: %f\n", inst->best_value);
            two_opt_refining(inst);
            printf("Two-opt cost cost: %f\n", inst->best_value);
            break;
        case XTRA_MILEAGE:
            printf("Model type chosen: undirected complete graph solved with extra mileage method\n");
            inst->start_time = seconds();
            extra_mileage(inst);
            printf("Solution cost: %f\n", inst->best_value);
            break;
        case XTRA_MILEAGE_REF:
            printf("Model type chosen: undirected complete graph solved with extra mileage method + "
                   "2-opt refining\n");
            inst->start_time = seconds();
            extra_mileage(inst);
            printf("Extra mileage cost: %f\n", inst->best_value);
            two_opt_refining(inst);
            printf("Two-opt cost cost: %f\n", inst->best_value);
            break;
        case VNS:
            printf("Model type chosen: VNS\n");
            inst->start_time = seconds();
            greedy(inst);
            vns(inst);
            printf("Solution cost: %f\n", inst->best_value);
            print_stats(inst, seconds()-inst->start_time);
            break;
        case TABU_SEARCH:
            printf("Model type chosen: Tabu Search\n");
            inst->start_time = seconds();
            greedy(inst);
            tabu_search(inst);
            printf("Solution cost: %f\n", inst->best_value);
            print_stats(inst, seconds()-inst->start_time);
            break;
        case GENETIC:
            printf("Model type chosen: Genetic Algorithm\n");
            inst->start_time = seconds();
            genetic(inst, 100);
            printf("Solution cost: %f\n", inst->best_value);
            print_stats(inst, seconds()-inst->start_time);
            break;
        default:
            TSPopt(inst);
            break;

    }
    double end_time = seconds();
    printf("Time elapsed: %f\n", end_time - start_time);
    return 0;
}

void performance_profile(instance *inst, int *models, int nmodels, double time_limit) {
    printf("------------START PERFORMANCE PROFILE MODE-----------\n");

    //Setting parameters inside of inst
    inst->timelimit = time_limit;
    printf("TIME LIMIT SETTED TO: %6.2fs\n", inst->timelimit);
    double start_time = seconds();

    FILE *csv = fopen("../logfiles/performance_profile.csv", "a");
    if (csv == NULL) {
        print_error_code("Failed to open ../logfiles/performance_profile.csv", ENOENT);
    }
    fprintf(csv, "%d, ", nmodels);

    //Print the first line of the csv file
    printf("Models selected:\n");
    for (int i = 0; i < nmodels; i++) {
        printf("%d\n", models[i]);
        fprintf(csv, "%d", models[i]);
        if (i != nmodels - 1)
            fprintf(csv, ", ");
        else
            fprintf(csv, "\n");
    }
    fclose(csv);

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
        csv = fopen("../logfiles/performance_profile.csv", "a");
        if (csv == NULL) {
            print_error_code("Failed to open ../logfiles/performance_profile.csv", ENOENT);
        } else {
            fprintf(csv, "%s, ", file_name);
            fclose(csv);
        }

        //Read and parse the file to test
        strcpy(inst->input_file, file_name);
        read_input(inst);

        //Iterate over the models passed
        for (int i = 0; i < nmodels; i++) {
            printf("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n");
            inst->model_type = models[i];
            double start_time_for = seconds();
            solve(inst);
            double end_time_for = seconds();

            csv = fopen("../logfiles/performance_profile_META.csv", "a");
            if (csv == NULL) {
                print_error_code("Failed to open ../logfiles/performance_profile.csv", ENOENT);
            }
            //Print
            if (inst->model_type == STANDARD || inst->model_type == BENDERS || inst->model_type == BRANCH_AND_CUT
                || inst->model_type == DEFAULT || inst->model_type == BRANCH_AND_CUT_RLX){
                fprintf(csv, "%f", end_time_for - start_time_for);
            } else {
                fprintf(csv, "%f", inst->best_value);
            }

            if (i != nmodels - 1)
                fprintf(csv, ", ");
            else
                fprintf(csv, "\n");
            fclose(csv);
            printf("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n\n");
            free_instance(inst, PERFORMANCE_PROFILE);
        }
        free(inst->latitude);
        free(inst->longitude);
        free(inst->x_coord);
        free(inst->y_coord);
        free(inst->component);
        free(inst->successors);
    }
    double end_time = seconds();
    printf("TIME ELAPSED: %f s\n\n", end_time - start_time);
}