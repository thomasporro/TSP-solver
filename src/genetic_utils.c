#include <tsp.h>
#include <utils.h>
#include <float.h>
#include "genetic_utils.h"

double generate_random_solution(instance *inst, int *node_list) {

    int *remaining_nodes = (int *) calloc(inst->nnodes, sizeof(int));
    for (int i = 0; i < inst->nnodes; i++) {
        remaining_nodes[i] = i;
    }

    int remaining_nodes_size = inst->nnodes - 1;
    int node_in_index = 0;

    while (remaining_nodes_size > 0) {
        int node_to_pop = rand() % remaining_nodes_size;
        node_list[node_in_index] = remaining_nodes[node_to_pop];
        remaining_nodes[node_to_pop] = remaining_nodes[remaining_nodes_size];
        node_in_index++;
        remaining_nodes_size--;
    }

    node_list[node_in_index] = remaining_nodes[remaining_nodes_size];

    double current_cost = distance(node_list[0], node_list[inst->nnodes - 1], inst);
    for (int i = 0; i < inst->nnodes - 1; i++) {
        current_cost += distance(node_list[i], node_list[i + 1], inst);
    }
    inst->best_value = current_cost;

    return current_cost;
}

void list_to_successors(instance *inst, const int *node_list, int *successors, const int *size) {
    int upper_limit;
    if (size != NULL) {
        upper_limit = *size;
    } else {
        upper_limit = inst->nnodes;
    }
    for (int i = 0; i < inst->nnodes; i++) {
        successors[i] = -1;
    }

    // Last to first
    successors[node_list[upper_limit - 1]] = node_list[0];
    for (int i = 0; i < upper_limit - 1; i++) {
        successors[node_list[i]] = node_list[i + 1];
    }
}

void successors_to_list(instance *inst, const int *successors, int *node_list) {
    int index = 0;
    int start_node = 0;
    int current_node = start_node;
    do {
        node_list[index] = successors[current_node];
        current_node = successors[current_node];
        index++;
    } while (current_node != start_node);
}

void complete_merging(instance *inst, int *node_list, int size) {
    int node_added = size;
    int start_node = node_list[0];
    int *successors = (int *) calloc(inst->nnodes, sizeof(int));
    list_to_successors(inst, node_list, successors, &size);

    while (node_added < inst->nnodes) {

        int current_node = start_node;
        int candidate_node = -1;
        int candidate_start = -1;
        int candidate_end = -1;
        double min_distance = DBL_MAX;

        //Find the node which extra mileage is the smallest among the free nodes
        do {
            int next_node = successors[current_node];

            for (int i = 0; i < inst->nnodes; i++) {
                if (successors[i] != -1 || i == current_node || i == next_node) continue;
                double delta = distance(current_node, i, inst) +
                               distance(i, next_node, inst) -
                               distance(current_node, next_node, inst);
                if (delta < min_distance) {
                    candidate_node = i;
                    candidate_start = current_node;
                    candidate_end = next_node;
                    min_distance = delta;
                }
            }

            current_node = next_node;
        } while (current_node != start_node);

        //Update the solution
        successors[candidate_start] = candidate_node;
        successors[candidate_node] = candidate_end;

        node_added++;
//        printf("node added: %d\n", node_added);
    }

    successors_to_list(inst, successors, node_list);
    free(successors);
}

double generate_greedy_solution(instance *inst, int *node_list) {
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    for (int i = 0; i < inst->nnodes; i++) {
        succ[i] = -1;
    }
    int start_node = rand() % inst->nnodes;

    int stop_flag = 0;

    int current_node = start_node;
    double cost = 0;

    while (!stop_flag) {
        double min_distance = DBL_MAX;
        int candidate_node = -1;

        for (int i = 0; i < inst->nnodes; i++) {
            //If the node analyzed is the current ore the successors of i has been already defined
            //skips the iteration
            if (current_node == i || succ[i] != -1) continue;

            //If the node of the cycle is nearest in respect the previous one change this values
            if (distance(current_node, i, inst) < min_distance) {
                min_distance = distance(current_node, i, inst);
                candidate_node = i;
            }
        }

        if (candidate_node != -1) {
            succ[current_node] = candidate_node;
            current_node = candidate_node;
            cost += min_distance;
        } else {
            succ[current_node] = start_node;
            cost += distance(current_node, start_node, inst);
            stop_flag = 1;
        }
    }

    successors_to_list(inst, succ, node_list);
    free(succ);
    return cost;
}

double generate_extra_mileage_solution(instance *inst, int *node_list) {
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    double best_value = -1;
    for (int i = 0; i < inst->nnodes; i++) {
        succ[i] = -1;
    }

    int node1 = rand() % inst->nnodes;
    int node2 = -1;
    while ((node2 = rand() % inst->nnodes) == node1);

    succ[node1] = node2;
    succ[node2] = node1;
    best_value = 2 * distance(node1, node2, inst);

    int start_node = node1;
    int node_added = 2;

    while (node_added < inst->nnodes) {

        int current_node = start_node;
        int candidate_node = -1;
        int candidate_start = -1;
        int candidate_end = -1;
        double min_distance = DBL_MAX;

        //Find the node which extra mileage is the smallest among the free nodes
        do {
            int next_node = succ[current_node];

            for (int i = 0; i < inst->nnodes; i++) {
                if (i == current_node || i == next_node || succ[i] != -1) continue;
                double delta = distance(current_node, i, inst) +
                               distance(i, next_node, inst) -
                               distance(current_node, next_node, inst);
                if (delta < min_distance) {
                    candidate_node = i;
                    candidate_start = current_node;
                    candidate_end = next_node;
                    min_distance = delta;
                }
            }

            current_node = next_node;
        } while (current_node != start_node);

        //Update the solution
        succ[candidate_start] = candidate_node;
        succ[candidate_node] = candidate_end;
        best_value = best_value -
                     distance(candidate_start, candidate_end, inst) +
                     distance(candidate_start, candidate_node, inst) +
                     distance(candidate_node, candidate_end, inst);

        node_added++;
    }
    successors_to_list(inst, succ, node_list);
    free(succ);
    return best_value;
}