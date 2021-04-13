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
	char edge_type[10];
	int model_type;

	double *latitude;
	double *longitude;

	//Variable that will contain global data
	double timelimit;
	double start_time;
	double *solution;
	int nvariables;
	int *successors;
	int *component;
	int ncomp;

} instance;

#endif // !TSP_H_


/*!
* Calculate the solution of the problem built into an instance
* @param	inst is a pointer to the instance where is stored the problem
* @return	0 if the solution is found. Other values otherwise
*/
int TSPopt(instance *inst);

/*!
* Function that switches the problem's build from a model to another using the type save into the instance
* @param	inst is a pointer to the instance of the problem created using tsp.h
* @param	env is the enviroment of CPLEX
* @param	lp is the problem written in CPLEX
*/
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);

/*!
* Function that build the model for an undirected graph
* @param	inst is a pointer to the instance of the problem created using tsp.h
* @param	env is the enviroment of CPLEX
* @param	lp is the problem written in CPLEX
*/
void build_model_st(instance *inst, CPXENVptr env, CPXLPptr lp);

/*!
* Function that build the MTZ model for an directed graph
* @param	inst is a pointer to the instance of the problem created using tsp.h
* @param	env is the enviroment of CPLEX
* @param	lp is the problem written in CPLEX
*/
void build_model_MTZ(instance *inst, CPXENVptr env, CPXLPptr lp);

/*!
* Function that build the GG model for an directed graph
* @param	inst is a pointer to the instance of the problem created using tsp.h
* @param	env is the enviroment of CPLEX
* @param	lp is the problem written in CPLEX
*/
void build_model_GG(instance *inst, CPXENVptr env, CPXLPptr lp);

/*!
* Function that build the solution for an undirected graph
* @param	inst is a pointer to the instance of the problem created using tsp.h
* @param	solution is the pointer to the solution
* @param	successors is the pointer to the array of successors
* @param	component is the pointer to the array where components are saved
* @param	ncomp is the pointer where to save the number of components
*/
void build_solution(instance *inst, double *solution, int *successors, int *component, int *ncomp);


/*!
* Function that start the CPXmipopt in many ways in order to execute 
* different methods
* @param	inst is a pointer to the instance of the problem created using tsp.h
* @param	env is the enviroment of CPLEX
* @param	lp is the problem written in CPLEX
*/
void compute_solution(instance *inst, CPXENVptr env, CPXLPptr lp);


/*!
* Implement the loop method by adding costraints for the loops found
* @param	inst is a pointer to the instance of the problem created using tsp.h
* @param	env is the enviroment of CPLEX
* @param	lp is the problem written in CPLEX
*/
void loop_method(instance *inst, CPXENVptr env, CPXLPptr lp);


/*!
* Cplex declartion of the user function to be called in the callback loop method
*/
static int CPXPUBLIC incumbent_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle);