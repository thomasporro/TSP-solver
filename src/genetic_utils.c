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
    int nnz=0;

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
    }

    successors_to_list(inst, successors, node_list);
}