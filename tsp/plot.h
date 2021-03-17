#include "tsp.h"

/*!
* Plot a graph using the passed command and the gnuplot enviroment
* @param	commands it's a vector of strings where are written the command of gnuplot
* @param	n_commands is an integer indicating the number of command passed into commands
* @param	inst is a pointer to the instance of the problem created using tsp.h
*/
void plot(char **commands, int n_commands, instance *inst);