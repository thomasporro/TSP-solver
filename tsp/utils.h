#include "tsp.h"


/*!
* Calculates the distance between two nodes using the euclidian distance
* @param	i is the index of the first node
* @param	j is the index of the second node
* @param	inst is a pointer to the instance where the nodes are stored
* @return	double that indicates the euclidian distance
*/
double distance(int i, int j, instance *inst);


/*!
* Calculates the position of the variable x(i, j) into the CPLEX problem for
* an undirected graph
* @param	i is the index of the first node
* @param	j is the index of the second node
* @param	inst is a pointer to the instance where the nodes are stored
* @return	an integer indicating the position of the variable into the CPLEX problem
*/
int xpos(int i, int j, instance *inst);


/*!
* Prints an error and exits from the program
* @param	err pointer to a char's string that will be print
*/
void print_error(const char *err);