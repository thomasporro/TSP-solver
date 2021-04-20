#define _CRT_SECURE_NO_DEPRECATE
#define BINARY_VARIABLE 'B'
#define INTEGER_VARIABLE 'I'
#define EQUAL 'E'
#define LESS_EQUAL 'L'
#define GREAT_EQUAL 'G'
#define EPS 1e-5

#include <stdint.h>
#include <time.h>
#include <concorde.h>
#include <cplex.h>
#include "../include/tsp.h"
#include "../include/utils.h"

int TSPopt(instance *inst) {
    //Open the CPLEX enviroment
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

    //Setting cplex parameters
    if (inst->timelimit != CPX_INFBOUND) {
        CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit);
    }

    CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);

    //Build the model into the cplex format
    build_model(inst, env, lp);

    //Obtain the number of variable and print them on screen
    int cur_numcols = CPXgetnumcols(env, lp);
    inst->nvariables = cur_numcols;
    printf("Number of variables: %d\n", inst->nvariables);

    inst->solution = (double *) calloc(inst->nvariables, sizeof(double));

    double start_time = seconds();

    //Computing the solution
    compute_solution(inst, env, lp);

    double time_passed = seconds() - start_time;
    print_stats(inst, time_passed);

    //If the problem have a solution it saves it into the instance's structure
    if (CPXgetx(env, lp, inst->solution, 0, inst->nvariables - 1)) {
        print_error("Failed to optimize MIP (retrieving the solution).\n");
    }

    //Print the values of the solution
    double solution = -1;
    if (CPXgetobjval(env, lp, &solution)) {
        print_error("Failed to optimize MIP.\n");
    }
    printf("objval: %f\n", solution);

    //Frees the memory
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return error;
}


void build_model(instance *inst, CPXENVptr env, CPXLPptr lp) {
    modeltype model = inst->model_type;
    switch (model) {
        case STANDARD:
            printf("Model type chosen: undirected complete graph\n");
            build_model_st(inst, env, lp);
            break;

        case BENDERS:
            printf("Model type chosen: undirected complete graph with loop method\n");
            build_model_st(inst, env, lp);
            break;

        case BRANCH_AND_CUT:
            printf("Model type chosen: undirected complete graph with loop method as a callback\n");
            build_model_st(inst, env, lp);
            break;

        case MTZ:
            printf("Model type chosen: MTZ\n");
            build_model_MTZ(inst, env, lp);
            break;

        case MTZ_LAZY:
            printf("Model type chosen: MTZ with lazy constraints\n");
            build_model_MTZ(inst, env, lp);
            break;

        case MTZ_IND:
            printf("Model type chosen: MTZ with indicator constraints\n");
            build_model_MTZ(inst, env, lp);
            break;

        case GG:
            printf("Model type chosen: GG\n");
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
    char **cname = (char **) calloc(1, sizeof(char *));
    cname[0] = (char *) calloc(100, sizeof(char));

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
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &variable_type, cname))
                print_error("Wrong CPXnewcols on x var.s");
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
    char **cname = (char **) calloc(1, sizeof(char *));
    cname[0] = (char *) calloc(100, sizeof(char));

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
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &variable_type, cname))
                print_error("Wrong CPXnewcols on x var.s");
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
        //TODO MANCA LA VARIABLE CHE FA IL CHECK
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
    if (inst->model_type == MTZ) { //Standard constraints
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
                if (CPXchgcoef(env, lp, lastrow, xxpos(i, j, inst), inst->nnodes - 1.0))
                    print_error("wrong CPXchgcoef [n*x(ij)]");
                if (CPXchgcoef(env, lp, lastrow, init_position_of_u + i, 1.0)) print_error("wrong CPXchgcoef [u(i)]");
                if (CPXchgcoef(env, lp, lastrow, init_position_of_u + j, -1.0)) print_error("wrong CPXchgcoef [u(j)]");
            }
        }
        // TODO: funziona anche se non � qui mi sa (CORRETTO!)
        /*
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
        }*/
    } else if (inst->model_type == MTZ_LAZY) {//Lazy constraints
        printf("Lazy constraints chosen correctly\n");
        int izero = 0;
        int index[3];
        double value[3];
        double big_m = inst->nnodes - 1.0;
        double rhs = big_m - 1.0;
        char sense = LESS_EQUAL;
        int nnz = 3;
        int final_position = inst->nnodes * inst->nnodes;
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
                if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname))
                    print_error("wrong lazy contraints for u-consistenxy");
            }
        }

        // add lazy constraints 1.0 * x_ij + 1.0 * x_ji <= 1, for each arc (i,j) with i < j
        /*
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
        }*/
    } else if (inst->model_type == MTZ_IND) { //Model solved with the indicators constraint
        printf("Indicator constraints chosen correctly\n");
        int rhs = 1;
        int nzcnt = 2;
        char sense = GREAT_EQUAL;
        int complemented = 0;
        int linind[2];
        double linval[] = {-1.0, 1.0};
        //Point the last position of the x_ij variables
        int final_position = inst->nnodes * inst->nnodes;
        for (int i = 1; i < inst->nnodes; i++) {
            for (int j = 1; j < inst->nnodes; j++) {
                if (i == j) continue;
                sprintf(cname[0], "indicator_constraint for arc(%d, %d)", i + 1, j + 1);
                linind[0] = final_position + i;
                linind[1] = final_position + j;
                if (CPXaddindconstr(env, lp, xxpos(i, j, inst), complemented, nzcnt, rhs, sense, linind, linval,
                                    cname[0])) {
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
    char **cname = (char **) calloc(1, sizeof(char *));
    cname[0] = (char *) calloc(100, sizeof(char));

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
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
            print_error("wrong CPXnewrows [flow_in_and_out_from_a_node]");
        for (int i = 0; i < inst->nnodes; i++) {
            if (j == i) continue;
            if (CPXchgcoef(env, lp, lastrow, ypos(i, j, inst), 1.0))
                print_error("wrong chgcoef [flow_in_and_out_from_a_node]");
            if (CPXchgcoef(env, lp, lastrow, ypos(j, i, inst), -1.0))
                print_error("wrong chgcoef [flow_in_and_out_from_a_node]");
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
        if (CPXchgcoef(env, lp, lastrow, xxpos(0, j, inst), -(inst->nnodes - 1)))
            print_error("wrong CPXchgcoef [linking_first_node]");
    }

    //For the other nodes
    for (int i = 1; i < inst->nnodes; i++) {
        for (int j = 0; j < inst->nnodes; j++) {
            if (i == j) continue;
            sprintf(cname[0], "linking_arc_x(%d,%d)", i + 1, j + 1);
            int lastrow = CPXgetnumrows(env, lp);
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("wrong CPXnewrows [linking_arc_x]");
            if (CPXchgcoef(env, lp, lastrow, ypos(i, j, inst), 1.0)) print_error("wrong chgcoef [linking_arc_x]");
            if (CPXchgcoef(env, lp, lastrow, xxpos(i, j, inst), -(inst->nnodes - 2)))
                print_error("wrong chgcoef [linking_arc_x]");
        }
    }

    printf("END OF BUILDING CONSTRAINTS\n");

    CPXwriteprob(env, lp, "model.lp", NULL);
    printf("END OF WRITING THE PROBLEM ON model.lp\n");
    free(cname[0]);
    free(cname);
}


void build_solution(instance *inst, double *solution, int *successors, int *component, int *ncomp) {
    //Verify if each node has exactly 2 edges
    int *node_degree = (int *) calloc(inst->nnodes, sizeof(int));
    //In this cycle I modify the degree of a node if the value of the edge is greater
    //than a tollerance value defined as EPS
    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = i + 1; j < inst->nnodes; j++) {

            //Verify the value of the edge between the nodes i and j
            int k = xpos(i, j, inst);
            if (fabs(solution[k]) > EPS && fabs(solution[k] - 1.0) > EPS) {
                print_error("wrong solution in build_solution()");
            }

            //Modify the values of the node's degree
            if (solution[k] > 0.5) {
                ++node_degree[i];
                ++node_degree[j];
            }
        }
    }

    //Check if the degree of each node is equal at 2
    for (int i = 0; i < inst->nnodes; i++) {
        if (node_degree[i] != 2) {
            printf("wrong degree in build_sol() in the node %d", i);
            exit(1);
        }
    }
    free(node_degree);

    //Let's set the values of successors and components and the number of components
    *ncomp = 0;
    for (int i = 0; i < inst->nnodes; i++) {
        successors[i] = -1;
        component[i] = -1;
    }

    //We build the successors and components to indicate the edges of the solution
    for (int start = 0; start < inst->nnodes; start++) {
        //If the component of start has been already assigned skip this for
        if (component[start] >= 0) continue;
        //I upgrade the number of components
        *ncomp = *ncomp + 1;
        int i = start;
        while (component[start] == -1) {
            for (int j = 0; j < inst->nnodes; j++) {
                if (i != j && solution[xpos(i, j, inst)] > 0.5
                    && component[j] == -1 && successors[j] != i) {
                    successors[i] = j;
                    component[j] = *ncomp;
                    i = j;
                    break;
                }
            }
        }
    }

    //printf("Number of loops: %d\n\n", *ncomp);
}


void compute_solution(instance *inst, CPXENVptr env, CPXLPptr lp) {
    modeltype model = inst->model_type;
    switch (model) {
        case BENDERS:
            inst->successors = (int *) calloc(inst->nnodes, sizeof(int));
            inst->component = (int *) calloc(inst->nnodes, sizeof(int));
            inst->ncomp = 99999;
            inst->start_time = seconds();

            int cycles = 1;
            while (inst->ncomp > 1 &&
                   seconds() - inst->start_time < inst->timelimit) {
                printf("CALCULATING THE SOLUTION (CYCLE N.%d) ...\n", cycles);
                CPXmipopt(env, lp);
                printf("SOLUTION CALCULATED (CYCLE N.%d)\n", cycles);
                printf("Number of costraints (CYCLE N.%d): %d\n", cycles, CPXgetnumrows(env, lp));

                //Change the internal time limit of cplex
                CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit - seconds());

                if (CPXgetx(env, lp, inst->solution, 0, inst->nvariables - 1)) {
                    print_error("Failed to optimize MIP.\n");
                }

                build_solution(inst, inst->solution, inst->successors, inst->component, &inst->ncomp);
                if (inst->ncomp >= 2) {
                    loop_method(inst, env, lp);
                }
                cycles++;
            }
            break;

        case BRANCH_AND_CUT:
            //Code needed to say to cplex that we have a callback function
            //TODO try to add the CPXLONG as a variable
            //CPXLONG contextId = CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION;
            if (CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION,
                                   sec_callback, inst)) {
                print_error("CPXcallbacksetfunc() error");
            }
            printf("CALCULATING THE SOLUTION...\n");
            CPXmipopt(env, lp);
            printf("SOLUTION CALCULATED\n");
            break;

        default:
            printf("CALCULATING THE SOLUTION...\n");
            CPXmipopt(env, lp);
            printf("SOLUTION CALCULATED\n");
            break;
    }
}


void loop_method(instance *inst, CPXENVptr env, CPXLPptr lp) {
    char **cname = (char **) calloc(1, sizeof(char *));
    cname[0] = (char *) calloc(100, sizeof(char));

    //Create a new costraint for each component
    for (int i = 1; i <= inst->ncomp; i++) {
        int *nodes = (int *) calloc(inst->nnodes, sizeof(int));
        int number_of_nodes = 0;
        for (int j = 0; j < inst->nnodes; j++) {
            if (inst->component[j] == i) {
                nodes[number_of_nodes] = j;
                number_of_nodes++;
            }
        }

        //Creates a new constraint and change the coef of all the nodes
        //that belongs to the current component
        int lastrow = CPXgetnumrows(env, lp);
        sprintf(cname[0], "SEC constraints");
        double rhs = number_of_nodes - 1;
        char sense = LESS_EQUAL;
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error(" wrong CPXnewrows [loop method]");

        for (int j = 0; j < number_of_nodes; j++) {
            for (int h = j + 1; h < number_of_nodes; h++) {
                if (j == h) continue;
                if (CPXchgcoef(env, lp, lastrow, xpos(nodes[j], nodes[h], inst), 1.0))
                    print_error("wrong CPXchgcoef [loop method]");
            }
        }
        free(nodes);
    }
    CPXwriteprob(env, lp, "model.lp", NULL);
    free(cname[0]);
    free(cname);
}


static int CPXPUBLIC sec_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle) {
    instance *inst = (instance *) userhandle;

    if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE) {
        return candidate_callback(context, inst);
    } else if (contextid == CPX_CALLBACKCONTEXT_RELAXATION) {
        return relaxation_callback(context, inst);
    }

    return 0;
}


int candidate_callback(CPXCALLBACKCONTEXTptr context, instance *inst) {
    //Retrieve the candidate solution
    double *solution_found = (double *) calloc(inst->nvariables, sizeof(double));
    double obj_value = CPX_INFBOUND;
    if (CPXcallbackgetcandidatepoint(context, solution_found, 0, inst->nvariables - 1, &obj_value)) {
        print_error("Error while retrieving the solution (incumbent_callback)");
    }


    int *successors = (int *) calloc(inst->nnodes, sizeof(int));
    int *component = (int *) calloc(inst->nnodes, sizeof(int));
    int ncomp = 99999;

    build_solution(inst, solution_found, successors, component, &ncomp);

    if (ncomp > 1) {
        //Iterate over the nodes of the component
        for (int i = 1; i <= ncomp; i++) {
            int *nodes = (int *) calloc(inst->nnodes, sizeof(int));
            int nzcnt = 0;
            for (int j = 0; j < inst->nnodes; j++) {
                if (component[j] == i) {
                    nodes[nzcnt] = j;
                    nzcnt++;
                }
            }

            //TODO try to improve the memory used -> placeit outside the for and give it a
            //predefined dimension
            //Build the arguments to pass to the function
            int *rmatind = (int *) calloc((nzcnt * nzcnt), sizeof(int));
            double *rmatval = (double *) calloc((nzcnt * nzcnt), sizeof(double));
            char sense = LESS_EQUAL;
            double rhs = nzcnt - 1;
            int rmatbeg = 0;

            int position = 0;
            for (int j = 0; j < nzcnt; j++) {
                for (int k = j + 1; k < nzcnt; k++) {
                    if (j == k) continue;
                    rmatind[position] = xpos(nodes[j], nodes[k], inst);
                    rmatval[position++] = 1.0;
                }
            }

            if (CPXcallbackrejectcandidate(context, 1, position, &rhs, &sense, &rmatbeg, rmatind, rmatval)) {
                print_error("Error on callbackrejectcandidate");
            }

            free(rmatind);
            free(rmatval);
            free(nodes);
        }
    }
    return 0;
}


int relaxation_callback(CPXCALLBACKCONTEXTptr context, instance *inst) {
    //Retrieve the fractional solution
    double *solution_found = (double *) calloc(inst->nvariables, sizeof(double));
    double obj_value = CPX_INFBOUND;
    if (CPXcallbackgetrelaxationpoint(context, solution_found, 0, inst->nvariables - 1, &obj_value)) {
        print_error("Error while retrieving the solution (incumbent_callback)");
    }

    //Inizialization of values for the concorde method
    int ncomp = 0;
    int *compscount = (int *) calloc(inst->nnodes, sizeof(int));
    int *comps = (int *) calloc(inst->nnodes, sizeof(int));
    int ecount = inst->nnodes * (inst->nnodes - 1) / 2;
    int *elist = (int *) calloc(ecount * 2, sizeof(int));
    int loader = 0;

    //Filling the elist array
    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = i + 1; j < inst->nnodes; j++) {
            //For each edge are reserved 2 slots over this array. Concorde
            //knows that the first edge are the first 2 slots and so on.
            //Now we're filling the array with the nodes of each edge
            elist[loader++] = i;
            elist[loader++] = j;
        }
    }

    if (CCcut_connect_components(inst->nnodes, ecount, elist, solution_found, &ncomp, &compscount, &comps)) {
        print_error("Error on CCcut_connect_components()");
    }

    if (ncomp != 1) {
        ccr_param inParam;
        inParam.context = context;
        inParam.inst = inst;
        if (CCcut_violated_cuts(inst->nnodes, ecount, elist, solution_found, 2.0 - EPS, create_cut_relaxation,
                                (void *)&inParam)) {
            print_error("Error on CCcut_violated_cuts()");
        }
    }
    return 0;
}


int create_cut_relaxation(double cutval, int cutcount, int *cut, void *inParam){
    ccr_param *parameters = (ccr_param *) inParam;

    //Build the arguments to pass to the function
    int *rmatind = (int *) calloc((cutcount * cutcount), sizeof(int));
    double *rmatval = (double *) calloc((cutcount * cutcount), sizeof(double));
    char sense = LESS_EQUAL;
    double rhs = cutcount - 1;
    int rmatbeg = 0;
    int purgeable = CPX_USECUT_FILTER;
    int local = 0;

    //The same as candidate_callback but the nodes are stored in *cut
    int position = 0;
    for (int j = 0; j < cutcount; j++) {
        for (int k = j + 1; k < cutcount; k++) {
            if (j == k) continue;
            rmatind[position] = xpos(cut[j], cut[k], parameters->inst);
            rmatval[position++] = 1.0;
        }
    }

    if(CPXcallbackaddusercuts(parameters->context, 1, cutcount, &rhs, &sense, &rmatbeg, rmatind, rmatval,
                              &purgeable, &local)){
        print_error("Error on CPXcallbackaddusercuts()");
    }

    return 0;
}