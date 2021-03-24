#include "tsp.h"


/*!
* Calculates the distance between two nodes using the euclidian distance
* @param	i is the index of the first node
* @param	j is the index of the second node
* @param	inst is a pointer to the instance where the nodes are stored
* @return	double that indicates the euclidian distance
*/
double euc_2d_distance(int i, int j, instance *inst);


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


/*!
* Calculates the position of the variable x(i, j) into the CPLEX problem for
* the MTZ and GG model
* @param	i is the index of the first node
* @param	j is the index of the second node
* @param	inst is a pointer to the instance where the nodes are stored
* @return	an integer indicating the position of the variable into the CPLEX problem
*/
int xxpos(int i, int j, instance *inst);

/*!
* Calculates the position of the variable y(i, j) into the CPLEX problem for
* GG model
* @param	i is the index of the first node
* @param	j is the index of the second node
* @param	inst is a pointer to the instance where the nodes are stored
* @return	an integer indicating the position of the variable into the CPLEX problem
*/
int ypos(int i, int j, instance *inst);


/*!
* Save into a file the time passed to solve an instance with a specific
* model type
* @param	inst is the pointer to the instance where we can retrieve informations
* @param	time_passed is the time passed to compute the solution in seconds
*/
void print_stats(instance *inst, double time_passed);


/*!
* Calculates the distance between two nodes using the geographical distance using
* the Haversine formula
* @param	i is the index of the first node
* @param	j is the index of the second node
* @param	inst is a pointer to the instance where the nodes are stored
* @return	double that indicates the euclidian distance
*/
double geo_distance(int i, int j, instance *inst);

/*!
* Calculates the distance between two nodes using the ATT type
* @param	i is the index of the first node
* @param	j is the index of the second node
* @param	inst is a pointer to the instance where the nodes are stored
* @return	double that indicates the euclidian distance
*/
double att_distance(int i, int j, instance *inst);


/*!
* Switch statement for the distance functions
* @param	i is the index of the first node
* @param	j is the index of the second node
* @param	inst is a pointer to the instance where the nodes are stored
* @return	double that indicates the euclidian distance
*/
double distance(int i, int j, instance *inst);
