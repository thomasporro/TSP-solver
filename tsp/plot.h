#include "tsp.h"

/*!
* Plot a graph using the passed command and the gnuplot enviroment
* @param	commands it's a vector of strings where are written the command of gnuplot
* @param	n_commands is an integer indicating the number of command passed into commands
* @param	inst is a pointer to the instance of the problem created using tsp.h
*/
void plot(char **commands, int n_commands, instance *inst);

/*!
* Print into the given file the solution of the undirected complete graph
* @param	temp It's the pointer to the file where we write the graph
* @param	inst is a pointer to the instance of the problem created using tsp.h
*/
void print_st(FILE *temp, instance *inst);

/*!
* Print into the given file the solution of the MTZ
* @param	temp It's the pointer to the file where we write the graph
* @param	inst is a pointer to the instance of the problem created using tsp.h
*/
void print_MTZ(FILE *temp, instance *inst);