#define _CRT_SECURE_NO_DEPRECATE
#define BINARY_VARIABLE 'B'
#define INTEGER_VARIABLE 'I'
#define EQUAL 'E'
#define LESS_EQUAL 'L'
#define GREAT_EQUAL 'G'

#include <stdint.h>
#include <time.h>
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
	
	time_t start_time = time(NULL);
	
	//Computing the solution
	printf("CALCULATING THE SOLUTION...\n");
	CPXmipopt(env, lp);
	printf("SOLUTION CALCULATED\n");

	double time_passed = difftime(time(NULL), start_time);
	print_stats(inst, time_passed);

	//If the problem have a solution it saves it into the instance's structure
	if (CPXgetx(env, lp, inst->solution, 0, inst->nvariables - 1)) {
		print_error("Failed to optimize MIP.\n");
	}

	//If the problem have a solution it saves it into the instance's structure
	double solution = -1;
	int lpstat = -1;
	if (CPXgetobjval(env, lp, &solution)){
		print_error("Failed to optimize MIP.\n");
	}
	printf("objval: %f\n", solution);

	//Frees the memory
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return error;
}


void build_model(instance *inst, CPXENVptr env, CPXLPptr lp) {
	switch (inst->model_type)
	{
	case 0:
		printf("Model type chosen: undirected complete graph\n");
		build_model_st(inst, env, lp);
		break;

	case 10:
		printf("Model type chosen: MTZ\n");
		build_model_MTZ(inst, env, lp);
		break;

	case 11:
		printf("Model type chosen: MTZ with lazy constraints\n");
		build_model_MTZ(inst, env, lp);
		break;

	case 12:
		printf("Model type chosen: MTZ with indicator constraints\n");
		build_model_MTZ(inst, env, lp);
		break;

	case 20:printf("Model type chosen: GG\n");
		build_model_GG(inst, env, lp);
		break;

	default:
		printf("Model type not defined. Will we used the model "
				"for the undirected complete graph\n");
		build_model_st(inst, env, lp);
		break;
	}
}


void build_model_st(instance *inst, CPXENVptr env, CPXLPptr lp) {
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


void build_model_MTZ(instance *inst, CPXENVptr env, CPXLPptr lp) {
	//Variables used to name the variables and constraints into CPLEX
	char **cname = (char **)calloc(1, sizeof(char*));
	cname[0] = (char *)calloc(100, sizeof(char));

	//Add binary variables x(i, j) for i < j
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = 0; j < inst->nnodes; j++) {
			//Variables used to create the variable into CPLEX
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
			char variable_type = BINARY_VARIABLE;
			double obj = distance(i, j, inst);
			double lb = 0.0;
			double ub = (i == j) ? 0.0 : 1.0;

			//Creates the columns of the new variable
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &variable_type, cname)) print_error("Wrong CPXnewcols on x var.s");
			if (CPXgetnumcols(env, lp) - 1 != xxpos(i, j, inst)) print_error("Wrong position for x var.s");
		}
	}

	//Add u(i) variable for MTZ
	for (int i = 0; i < inst->nnodes; i++) {
		sprintf(cname[0], "u(%d)", i + 1);
		char integer = INTEGER_VARIABLE;
		double lb = 0.0;
		double ub = inst->nnodes - 2.0;
		if (CPXnewcols(env, lp, 1, NULL, &lb, &ub, &integer, cname)) print_error("wrong CPXnewcols on x var.s");
	}

	//Add degree IN
	for (int h = 0; h < inst->nnodes; h++) {
		int lastrow = CPXgetnumrows(env, lp);
		double rhs = 1.0;
		char sense = EQUAL;
		sprintf(cname[0], ("degree_in_node_(%d)"), h + 1);
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(" wrong CPXnewrows [degree IN]");
		for (int j = 0; j < inst->nnodes; j++) {
			if (CPXchgcoef(env, lp, lastrow, xxpos(j, h, inst), 1.0)) print_error("wrong CPXchgcoef [degree IN]");
		}
	}

	//Add degree OUT
	for (int h = 0; h < inst->nnodes; h++) {
		int lastrow = CPXgetnumrows(env, lp);
		double rhs = 1.0;
		char sense = EQUAL;
		sprintf(cname[0], ("degree_out_node_(%d)"), h + 1);
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(" wrong CPXnewrows [degree OUT]");
		for (int j = 0; j < inst->nnodes; j++) {
			if (CPXchgcoef(env, lp, lastrow, xxpos(h, j, inst), 1.0)) print_error("wrong CPXchgcoef [degree OUT]");
		}
	}

	//Choose which constraint it should use
	if (inst->model_type == 10) { //Standard constraints
		printf("Standard constraints chosen correctly\n");
		//I skipped the first node since the model requires it
		for (int i = 1; i < inst->nnodes; i++) {
			for (int j = 1; j < inst->nnodes; j++) {
				if (i == j) continue;
				int lastrow = CPXgetnumrows(env, lp);
				double rhs = inst->nnodes - 2.0;
				char sense = LESS_EQUAL;

				//Compute the position of the first variable ui, that is right after the binary ones
				int init_position_of_u = xxpos(inst->nnodes - 1, inst->nnodes - 1, inst) + 1;
				sprintf(cname[0], ("MTZ_constraint_for_x(%d,%d)"), i + 1, j + 1);
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("wrong CPXnewrows [MTZ Constraint]");
				if (CPXchgcoef(env, lp, lastrow, xxpos(i, j, inst), inst->nnodes - 1.0)) print_error("wrong CPXchgcoef [n*x(ij)]");
				if (CPXchgcoef(env, lp, lastrow, init_position_of_u + i, 1.0)) print_error("wrong CPXchgcoef [u(i)]");
				if (CPXchgcoef(env, lp, lastrow, init_position_of_u + j, -1.0)) print_error("wrong CPXchgcoef [u(j)]");
			}
		}
		for (int i = 1; i < inst->nnodes; i++) {
			for (int j = i + 1; j < inst->nnodes; j++) {
				if (i == j) continue;
				int lastrow = CPXgetnumrows(env, lp);
				double rhs = 1.0;
				char sense = LESS_EQUAL;

				sprintf(cname[0], ("MTZ_constraint_for_xij_xji_1(%d,%d)"), i + 1, j + 1);
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("wrong CPXnewrows [MTZ Constraint]");
				if (CPXchgcoef(env, lp, lastrow, xxpos(i, j, inst), 1.0)) print_error("wrong CPXchgcoef [xij+xji<=1]");
				if (CPXchgcoef(env, lp, lastrow, xxpos(j, i, inst), 1.0)) print_error("wrong CPXchgcoef [xij+xji<=1]");
			}
		}
	}
	else if(inst->model_type == 11){//Lazy constraints
		printf("Lazy constraints chosen correctly\n");
		int izero = 0;
		int index[3];
		double value[3];
		double big_m = inst->nnodes - 1.0;
		double rhs = big_m - 1.0;
		char sense = LESS_EQUAL;
		int nnz = 3;
		int final_position = inst->nnodes  * inst->nnodes;
		for (int i = 1; i < inst->nnodes; i++) {
			for (int j = 1; j < inst->nnodes; j++) {
				if (i == j) continue;
				sprintf(cname[0], "u_consistency for arc (%d,%d)", i + 1, j + 1);
				index[0] = xxpos(i, j, inst);
				index[1] = final_position + i;
				index[2] = final_position + j;
				value[0] = big_m;
				value[1] = 1.0;
				value[2] = -1.0;
				if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname)) print_error("wrong lazy contraints for u-consistenxy");
			}
		}

		// add lazy constraints 1.0 * x_ij + 1.0 * x_ji <= 1, for each arc (i,j) with i < j
		rhs = 1.0;
		sense = LESS_EQUAL;
		nnz = 2;
		for (int i = 0; i < inst->nnodes; i++)
		{
			for (int j = i + 1; j < inst->nnodes; j++)
			{
				sprintf(cname[0], "SEC on node pair (%d,%d)", i + 1, j + 1);
				index[0] = xxpos(i, j, inst);
				value[0] = 1.0;
				index[1] = xxpos(j, i, inst);
				value[1] = 1.0;
				if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname)) print_error("wrong CPXlazyconstraints on 2-node SECs");
			}
		}
	}
	else if (inst->model_type == 12) { //Model solved with the indicators constraint
		printf("Indicator constraints chosen correctly\n");
		int rhs = 1;
		int nzcnt = 2;
		char sense = GREAT_EQUAL;
		int complemented = 0;
		int linind[2];
		double linval[] = { -1.0, 1.0};
		//Point the last position of the x_ij variables
		int final_position = inst->nnodes  * inst->nnodes;
		for (int i = 1; i < inst->nnodes; i++) {
			for (int j = 1; j < inst->nnodes; j++) {
				if (i == j) continue;
				sprintf(cname[0], "indicator_constraint for arc(%d, %d)", i+1, j+1);
				linind[0] = final_position + i;
				linind[1] = final_position + j;
				if (CPXaddindconstr(env, lp, xxpos(i, j, inst), complemented, nzcnt, rhs, sense, linind, linval, cname[0])) {
					print_error("wrong indicator contraint");
				}
			}
		}
	}
	


	printf("END OF BUILDING CONSTRAINTS\n");

	CPXwriteprob(env, lp, "model.lp", NULL);
	printf("END OF WRITING THE PROBLEM ON model.lp\n");
	free(cname[0]);
	free(cname);
}


void build_model_GG(instance *inst, CPXENVptr env, CPXLPptr lp) {
	char **cname = (char **)calloc(1, sizeof(char*));
	cname[0] = (char *)calloc(100, sizeof(char));

	//Add binary var.s x(i, j)
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = 0; j < inst->nnodes; j++) {
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
			char binary = BINARY_VARIABLE;
			double obj = distance(i, j, inst);
			double lb = 0.0;
			double ub = (i == j) ? 0.0 : 1.0;
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) print_error("wrong CPXnewcols on x var.s");
			if (CPXgetnumcols(env, lp) - 1 != xxpos(i, j, inst)) print_error("wrong position for x var.s");
		}
	}

	//TODO undestand ypos from the last year
	//Add variable y(i, j) corresponding to the flow of the arc
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = 0; j < inst->nnodes; j++) {
			sprintf(cname[0], "y(%d,%d)", i + 1, j + 1);
			char integer = INTEGER_VARIABLE;
			double lb = 0.0;
			double ub = (i == j) ? 0.0 : INFINITY;
			if (CPXnewcols(env, lp, 1, NULL, &lb, &ub, &integer, cname)) print_error("wrong CPXnewcols on y(i, j)");
			if (CPXgetnumcols(env, lp) - 1 != ypos(i, j, inst)) print_error("wrong position for y(i, j)");
		}
	}

	//Add degree IN
	for (int h = 0; h < inst->nnodes; h++) {
		int lastrow = CPXgetnumrows(env, lp);
		double rhs = 1.0;
		char sense = EQUAL;
		sprintf(cname[0], ("degree_in_node_(%d)"), h + 1);
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("wrong CPXnewrows [degree IN]");
		for (int j = 0; j < inst->nnodes; j++) {
			if (CPXchgcoef(env, lp, lastrow, xxpos(j, h, inst), 1.0)) print_error("wrong CPXchgcoef [degree IN]");
		}
	}

	//Add degree OUT
	for (int h = 0; h < inst->nnodes; h++) {
		int lastrow = CPXgetnumrows(env, lp); 
		double rhs = 1.0;
		char sense = EQUAL;
		sprintf(cname[0], ("degree_out_node_(%d)"), h + 1);
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("wrong CPXnewrows [degree OUT]");
		for (int j = 0; j < inst->nnodes; j++) {
			if (CPXchgcoef(env, lp, lastrow, xxpos(h, j, inst), 1.0)) print_error("wrong CPXchgcoef [degree OUT]");
		}
	}

	//Add costraint of flow of the first node
	sprintf(cname[0], "flow_out_from_1");
	int lastrow = CPXgetnumrows(env, lp);
	double rhs = inst->nnodes - 1;
	char sense = EQUAL;
	if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("wrong CPXnewrows [flow_out_from_1]");
	for (int j = 1; j < inst->nnodes; j++) {
		if (CPXchgcoef(env, lp, lastrow, ypos(0, j, inst), 1.0)) print_error("wrong chgcoef [flow_out_from_1]");
	}

	//Add costraint of flow of the other nodes
	for (int j = 1; j < inst->nnodes; j++) {
		sprintf(cname[0], "flow_in_and_out_from(%d)", j + 1);
		int lastrow = CPXgetnumrows(env, lp);
		double rhs = 1.0;
		char sense = EQUAL;
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("wrong CPXnewrows [flow_in_and_out_from_a_node]");
		for (int i = 0; i < inst->nnodes; i++) {
			if (j == i) continue;
			if (CPXchgcoef(env, lp, lastrow, ypos(i, j, inst), 1.0)) print_error("wrong chgcoef [flow_in_and_out_from_a_node]");
			if (CPXchgcoef(env, lp, lastrow, ypos(j, i, inst), -1.0)) print_error("wrong chgcoef [flow_in_and_out_from_a_node]");
		}
	}

	//Linking constraints 
	sprintf(cname[0], "linking_first_node");
	rhs = 0;
	sense = LESS_EQUAL;

	//For the first node
	for (int j = 1; j < inst->nnodes; j++) {
		int lastrow = CPXgetnumrows(env, lp);
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("wrong CPXnewrows [linking_first_node]");
		if (CPXchgcoef(env, lp, lastrow, ypos(0, j, inst), 1.0)) print_error("wrong CPXchgcoef [linking_first_node]");
		if (CPXchgcoef(env, lp, lastrow, xxpos(0, j, inst), -(inst->nnodes-1))) print_error("wrong CPXchgcoef [linking_first_node]");
	}

	//For the other nodes
	for (int i = 1; i < inst->nnodes; i++) {
		for (int j = 1; j < inst->nnodes; j++) {
			if (i == j) continue;
			sprintf(cname[0], "linking_arc_x(%d,%d)", i+1, j+1);
			int lastrow = CPXgetnumrows(env, lp); 
			if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("wrong CPXnewrows [linking_arc_x]");
			if (CPXchgcoef(env, lp, lastrow, ypos(i, j, inst), 1.0)) print_error("wrong chgcoef [linking_arc_x]");
			if (CPXchgcoef(env, lp, lastrow, xxpos(i, j, inst), -(inst->nnodes - 2))) print_error("wrong chgcoef [linking_arc_x]");
		}
	}

	printf("END OF BUILDING CONSTRAINTS\n");

	CPXwriteprob(env, lp, "model.lp", NULL);
	printf("END OF WRITING THE PROBLEM ON model.lp\n");
	free(cname[0]);
	free(cname);
}
