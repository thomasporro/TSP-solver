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
* Prints an error and exits from the program
* @param	err pointer to a char's string that will be print
* @param    status the code of the error
*/
void print_error_code(const char *err, int status);


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


/*!
* Return the seconds of the current program 
*/
double seconds();


/*!
 * Create and save the path name of the logfile for cplex
 * @param inst The instance, it is used in order to retrieve the model type
 * @return The array where the string will be written
 */
char *logfilename(instance *inst);


/*!
 * Frees the memory of the instance passed as argument
 * @param inst The instance that we want to kill
 * @param model_type 1 if we want to free the memory occupied by the solution of cplex
 */
void free_instance(instance *inst, int model_type);


/*!
 * Compute the solution starting from the values of successors
 * @param inst The instance
 * @param x The solution to be saved
 * @param successors The array containing the successors
 */
void compute_solution_from_successors(instance *inst, double *x, int *successors);


/*!
 * Compute the solution starting from the values of successors
 * @param inst The instance
 * @param successors The array containing the successors
 * @return The cost of the solution
 */
double compute_solution_cost(instance *inst, int *successors);


/*!
 * Compute the best 2-cut to perform over a solution
 * @param inst The instace where the solution is saved
 * @param first_node The first node found
 * @param second_node The second node found
 * @param improvement The improvement of the solution
 */
void compute_bigger_cut(instance *inst, int *first_node, int *second_node, double *improvement);


/*!
 * Performs the 2-opt refining move
 * @param inst The instance containing the solution
 * @param first_node The first node of the cut
 * @param second_node The second node of the cut
 * @param improvement The total improvement applied to the solution
 */
void perform_cut(instance *inst, int first_node, int second_node, double improvement);


