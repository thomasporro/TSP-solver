#define _CRT_SECURE_NO_DEPRECATE
#define BINARY_VARIABLE 'B'
#define INTEGER_VARIABLE 'I'
#define EQUAL 'E'
#define LESS_EQUAL 'L'
#define GREAT_EQUAL 'G'

#include <ilcplex/cplex.h>
#include "tsp.h"
#include "utils.h"

int TSPopt(instance *inst) {
	//Open the CPLEX enviroment
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

	//Build the model into the cplex format
	build_model(inst, env, lp);

	//Obtain the number of variable and print them on screen
	int cur_numcols = CPXgetnumcols(env, lp);
	inst->nvariables = cur_numcols;
	printf("Number of variables: %d\n", inst->nvariables);


	inst->solution = (double *)calloc(inst->nvariables, sizeof(double));

	
	//Computing the solution
	printf("CALCULATING THE SOLUTION...\n");
	CPXmipopt(env, lp);
	printf("SOLUTION CALCULATED\n");

	//If the problem have a solution it saves it into the instance's structure
	if (CPXgetx(env, lp, inst->solution, 0, inst->nvariables - 1)) {
		print_error("Failed to optimize MIP.\n");
	}


	//Frees the memory
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return error;
}


void build_model(instance *inst, CPXENVptr env, CPXLPptr lp) {
	//Variables used to name the variables and constraints into CPLEX
	char **cname = (char **)calloc(1, sizeof(char*));
	cname[0] = (char *)calloc(100, sizeof(char));

	//Add binary variables x(i, j) for i < j
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = i + 1; j < inst->nnodes; j++) {
			//Variables used to create the variable into CPLEX
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
			char variable_type = BINARY_VARIABLE;
			double obj = distance(i, j, inst);
			double lb = 0.0;
			double ub = 1.0;

			//Creates the columns of the new variable
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &variable_type, cname)) print_error("Wrong CPXnewcols on x var.s");
			if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, inst)) print_error("Wrong position for x var.s");
		}
	}

	//Add degree constraints for each node. The graph is undirected so the
	//degree for each node is 2.
	for (int h = 0; h < inst->nnodes; h++) {
		//Gives the row where to write the constraint
		int lastrow = CPXgetnumrows(env, lp);

		//Variables of the constraints
		double rhs = 2.0;
		char sense = EQUAL;
		sprintf(cname[0], "Degree(%d)", h + 1);

		//Adds the row of the constraint
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("Wrong CPXnewrows [degree]");

		//Modify the coefficients of the constrain putting the desidered value for the wanted variable
		for (int i = 0; i < inst->nnodes; i++) {
			if (i == h) continue;
			if (CPXchgcoef(env, lp, lastrow, xpos(i, h, inst), 1.0)) print_error("Wrong CPXchgcoef [degree]");
		}
	}


	printf("END OF BUILDING CONSTRAINTS\n");

	CPXwriteprob(env, lp, "model.lp", NULL);
	printf("END OF WRITING THE PROBLEM ON model.lp\n");
	free(cname[0]);
	free(cname);
}