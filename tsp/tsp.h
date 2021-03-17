#include <ilcplex/cplex.h>

#ifndef TSP_H_
#define TSP_H_

//Struct that will contain the parameters of the tsp problem
typedef struct {

	//Variables of the input files
	int nnodes;
	double *x_coord;
	double *y_coord;
	char input_file[1000];
	int model_type;

	//Variable that will contain global data
	double timelimit;
	double *solution;
	int nvariables;

} instance;

#endif // !TSP_H_


/*!
* Calculate the solution of the problem built into an instance
* @param	inst is a pointer to the instance where is stored the problem
* @return	0 if the solution is found. Other values otherwise
*/
int TSPopt(instance *inst);

/*!
* Function that build the model for an undirected graph
* @param	inst is a pointer to the instance of the problem created using tsp.h
* @param	env is the enviroment of CPLEX
* @param	lp is the problem written in CPLEX
*/
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);