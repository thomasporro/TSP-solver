#ifndef GENETIC_UTILS
#define GENETIC_UTILS

#include <stdio.h>
#include <string.h>

/*!
 * Generate a random feasible solution
 * @param inst The instance containing basic information
 * @param node_list A list of nodes that forms the solution. The nodes are in sequence
 * @return The cost of the random solution
 */
double generate_random_solution(instance *inst, int *node_list);

/*!
 * From a sequence of nodes generates the successors array
 * @param inst The instance containing basic information
 * @param node_list The sequence of nodes
 * @param successors Output array
 */
void list_to_successors(instance *inst, const int *node_list, int *successors, const int *size);

/*!
 * From the successors array generates a sequence of nodes
 * @param inst he instance containing basic information
 * @param successors The solution as successors
 * @param node_list Output array
 */
void successors_to_list(instance *inst, const int *successors, int *node_list);

/*!
 * Complete the merging of the chromosomes by using the extra-mileage algorithm
 * @param inst The instance containing basic information
 * @param node_list The list where the nodes are contained, it will be also the output array
 * @param size The size of the initial node_list array
 */
void complete_merging(instance *inst, int *node_list, int size);

double generate_greedy_solution(instance *inst, int *node_list);

double generate_extra_mileage_solution(instance *inst, int *node_list);

#endif //GENETIC_UTILS
