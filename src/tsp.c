#define BINARY_VARIABLE 'B'
#define INTEGER_VARIABLE 'I'
#define EQUAL 'E'
#define LESS_EQUAL 'L'
#define GREAT_EQUAL 'G'
#define LOWER_BOUND 'L'
#define EPS 1e-5
#define INTERNAL_TIME_LIMIT 60.0

#include <time.h>
#include <concorde.h>
#include <cplex.h>
#include <float.h>
#include "tsp.h"
#include "utils.h"
#include "genetic_utils.h"

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

    if (inst->performance_profile) {
        printf("-----Setting parameters for performance profile-----\n");
        //Setting the saving of the node's tree on the memory
        if (CPXsetstrparam(env, CPX_PARAM_WORKDIR, "../logfiles")) {
            print_error("CPX_PARAM_WORKDIR not setted");
        }
        if (CPXsetintparam(env, CPX_PARAM_NODEFILEIND, 2)) {
            print_error("CPX_PARAM_NODEFILEIND not setted");
        }
    }

    CPXsetlogfilename(env, logfilename(inst), "w");

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

    if (inst->model_type != HARD_FIX_BAC && inst->model_type != SOFT_FIX) {
        //If the problem have a solution it saves it into the instance's structure
        int status = CPXgetx(env, lp, inst->solution, 0, inst->nvariables - 1);
        if (status) {
            print_error_code("Failed to optimize MIP (retrieving the solution -> TSPopt())", status);
        }

        //Print the values of the solution
        double solution = -1;
        if (CPXgetobjval(env, lp, &solution)) {
            print_error("Failed to optimize MIP (getobjval() -> TSPopt).\n");
        }
        printf("Objval: %f\n", solution);
    }

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

        case BRANCH_AND_CUT_RLX:
            printf("Model type chosen: undirected complete graph with loop method as a callback "
                   "with relaxation\n");
            build_model_st(inst, env, lp);
            break;

        case HARD_FIX_BAC:
            printf("Model type chosen: callback + hard fixing\n");
            build_model_st(inst, env, lp);
            break;

        case SOFT_FIX:
            printf("Model type chosen: soft fixing\n");
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

    CPXwriteprob(env, lp, "../logfiles/model_st.lp", NULL);
    printf("END OF WRITING THE PROBLEM ON logfiles/model_st.lp\n");
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
    int start_position_u_i = xxpos(inst->nnodes - 1, inst->nnodes - 1, inst);
    for (int i = 1; i < inst->nnodes; i++) {
        sprintf(cname[0], "u(%d)", i + 1);
        char integer = INTEGER_VARIABLE;
        double lb = 0.0;
        double ub = inst->nnodes - 2.0;
        if (CPXnewcols(env, lp, 1, NULL, &lb, &ub, &integer, cname)) print_error("wrong CPXnewcols on x var.s");
        if (CPXgetnumcols(env, lp) - 1 != start_position_u_i + i) print_error("Wrong position for u");
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
                int init_position_of_u = xxpos(inst->nnodes - 1, inst->nnodes - 1, inst);
                sprintf(cname[0], ("MTZ_constraint_for_x(%d,%d)"), i + 1, j + 1);
                if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("wrong CPXnewrows [MTZ Constraint]");
                if (CPXchgcoef(env, lp, lastrow, xxpos(i, j, inst), inst->nnodes - 1.0))
                    print_error("wrong CPXchgcoef [n*x(ij)]");
                if (CPXchgcoef(env, lp, lastrow, init_position_of_u + i, 1.0)) print_error("wrong CPXchgcoef [u(i)]");
                if (CPXchgcoef(env, lp, lastrow, init_position_of_u + j, -1.0)) print_error("wrong CPXchgcoef [u(j)]");
            }
        }

    } else if (inst->model_type == MTZ_LAZY) {//Lazy constraints
        printf("Lazy constraints chosen correctly\n");
        int izero = 0;
        int index[3];
        double value[3];
        double big_m = inst->nnodes - 1.0;
        double rhs = big_m - 1.0;
        char sense = LESS_EQUAL;
        int nnz = 3;
        int final_position = xxpos(inst->nnodes - 1, inst->nnodes - 1, inst);
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

    } else if (inst->model_type == MTZ_IND) { //Model solved with the indicators constraint
        printf("Indicator constraints chosen correctly\n");
        int rhs = 1;
        int nzcnt = 2;
        char sense = GREAT_EQUAL;
        int complemented = 0;
        int linind[2];
        double linval[] = {-1.0, 1.0};
        //Point the last position of the x_ij variables
        int final_position = xxpos(inst->nnodes - 1, inst->nnodes - 1, inst);
        for (int i = 1; i < inst->nnodes; i++) {
            for (int j = 1; j < inst->nnodes; j++) {
                if (i == j) continue;
                sprintf(cname[0], "indicator_constraint for arc(%d, %d)", i + 1, j + 1);
                linind[0] = final_position + i;
                linind[1] = final_position + j;
                if (CPXaddindconstr(env, lp, xxpos(i, j, inst), complemented, nzcnt, rhs, sense, linind, linval,
                                    cname[0])) {
                    print_error("wrong indicator constraint");
                }
            }
        }
    }


    printf("END OF BUILDING CONSTRAINTS\n");

    CPXwriteprob(env, lp, "../logfiles/model_mtz.lp", NULL);
    printf("END OF WRITING THE PROBLEM ON logfiles/model_mtz.lp\n");
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
            double ub = (i == j) ? 0.0 : inst->nnodes - 1;
            if (CPXnewcols(env, lp, 1, NULL, &lb, &ub, &integer, cname)) print_error("wrong CPXnewcols on y(i, j)");
            if (CPXgetnumcols(env, lp) - 1 != ypos(i, j, inst)) print_error("wrong position for y(i, j)");
        }
    }

    //Variables used during the constraints
    int lastrow;
    double rhs;
    char sense;

    //Add degree IN
    for (int h = 0; h < inst->nnodes; h++) {
        lastrow = CPXgetnumrows(env, lp);
        rhs = 1.0;
        sense = EQUAL;
        sprintf(cname[0], ("degree_in_node_(%d)"), h + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("wrong CPXnewrows [degree IN]");
        for (int j = 0; j < inst->nnodes; j++) {
            if (CPXchgcoef(env, lp, lastrow, xxpos(j, h, inst), 1.0)) print_error("wrong CPXchgcoef [degree IN]");
        }
    }

    //Add degree OUT
    for (int h = 0; h < inst->nnodes; h++) {
        lastrow = CPXgetnumrows(env, lp);
        rhs = 1.0;
        sense = EQUAL;
        sprintf(cname[0], ("degree_out_node_(%d)"), h + 1);
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("wrong CPXnewrows [degree OUT]");
        for (int j = 0; j < inst->nnodes; j++) {
            if (CPXchgcoef(env, lp, lastrow, xxpos(h, j, inst), 1.0)) print_error("wrong CPXchgcoef [degree OUT]");
        }
    }

    //Add constraint of flow of the first node
    sprintf(cname[0], "flow_out_from_1");
    lastrow = CPXgetnumrows(env, lp);
    rhs = inst->nnodes - 1.0;
    sense = EQUAL;
    if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("wrong CPXnewrows [flow_out_from_1]");
    for (int j = 1; j < inst->nnodes; j++) {
        if (CPXchgcoef(env, lp, lastrow, ypos(0, j, inst), 1.0)) print_error("wrong chgcoef [flow_out_from_1]");
    }

    //Add constraint of flow of the other nodes
    for (int j = 1; j < inst->nnodes; j++) {
        sprintf(cname[0], "flow_in_and_out_from(%d)", j + 1);
        lastrow = CPXgetnumrows(env, lp);
        rhs = 1.0;
        sense = EQUAL;
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


    //Tightened linking constraints
    double coef_x = -(inst->nnodes - 2.0);
    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = 0; j < inst->nnodes; j++) {
            if (i == j) continue;
            if (i == 0 || j == 0) { coef_x = -(inst->nnodes - 1.0); }
            else { coef_x = -(inst->nnodes - 2.0); }
            sprintf(cname[0], "linking_arc_x(%d,%d)", i + 1, j + 1);
            lastrow = CPXgetnumrows(env, lp);
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) print_error("wrong CPXnewrows [linking_arc_x]");
            if (CPXchgcoef(env, lp, lastrow, ypos(i, j, inst), 1.0)) print_error("wrong chgcoef [linking_arc_x]");
            if (CPXchgcoef(env, lp, lastrow, xxpos(i, j, inst), coef_x))
                print_error("wrong chgcoef [linking_arc_x]");
        }
    }

    printf("END OF BUILDING CONSTRAINTS\n");

    CPXwriteprob(env, lp, "../logfiles/model_gg.lp", NULL);
    printf("END OF WRITING THE PROBLEM ON logfiles/model_gg.lp\n");
    free(cname[0]);
    free(cname);
}


void build_solution(instance *inst, double *solution, int *successors, int *component, int *ncomp) {
    //Verify if each node has exactly 2 edges
    int *node_degree = (int *) calloc(inst->nnodes, sizeof(int));
    //In this cycle I modify the degree of a node if the value of the edge is greater
    //than a tolerance value defined as EPS
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

                //If the time remaining is lower than the INTERNAL_TIME_LIMIT than
                //change the cplex time limit to the actual time remaining
                CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit - (seconds() - inst->start_time));


                printf("CALCULATING THE SOLUTION (CYCLE N.%d) ...\n", cycles);
                double cycle_time_start = seconds();
                CPXmipopt(env, lp);
                double cycle_time_end = seconds();
                printf("SOLUTION CALCULATED (CYCLE N.%d)\n", cycles);
                printf("Number of costraints (CYCLE N.%d): %d\n", cycles, CPXgetnumrows(env, lp));

                //Change the internal time limit of cplex
                CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit - seconds());

                if (CPXgetx(env, lp, inst->solution, 0, inst->nvariables - 1)) {
                    print_error("Failed to optimize MIP (getx() -> compute_solution()).\n");
                }

                build_solution(inst, inst->solution, inst->successors, inst->component, &inst->ncomp);
                printf("Number of loops: %d\n", inst->ncomp);
                if (inst->ncomp >= 2) {
                    loop_method(inst, env, lp);
                }
                cycles++;
                printf("Cycle time consumed: %f\n", cycle_time_end - cycle_time_start);
                printf("Total time consumed: %f\n", seconds() - inst->start_time);
                printf("---------------------------------------\n");
            }
            break;

        case BRANCH_AND_CUT:
            //Code needed to say to cplex that we have a callback function
            if (CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE, sec_callback, inst)) {
                print_error("CPXcallbacksetfunc() error -> BRANCH_AND_CUT");
            }
            printf("CALCULATING THE SOLUTION...\n");
            CPXmipopt(env, lp);
            printf("SOLUTION CALCULATED\n");
            break;

        case BRANCH_AND_CUT_RLX:
            //Code needed to say to cplex that we have a callback function
            if (CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION,
                                   sec_callback, inst)) {
                print_error("CPXcallbacksetfunc() error -> BRANCH_AND_CUT_RLX");
            }
            printf("CALCULATING THE SOLUTION...\n");
            CPXmipopt(env, lp);
            printf("SOLUTION CALCULATED\n");
            break;

        case HARD_FIX_BAC:
            inst->best_value = INFINITY;
            //Code needed to say to cplex that we have a callback function
            if (CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE,
                                   sec_callback, inst)) {
                print_error("CPXcallbacksetfunc() error -> HARD_FIX");
            }

            inst->start_time = seconds();
            //Change the internal time limit of cplex
            CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit / 3);

            //Setting up the array of indices for changing the bounds
            int *indices = (int *) calloc(inst->nnodes * 2, sizeof(int));
            char *lu = (char *) calloc(inst->nnodes * 2, sizeof(char));
            double *bd = (double *) calloc(inst->nnodes * 2, sizeof(double));

            //Counts the number of variables fixed in each cycle
            int cnt = 0;

            //Percentage with which I block the edges that are set to 1
            int percentage = 90;

            printf("First solving\n");
            CPXmipopt(env, lp);

            printf("Percentage of edges blocked: %d%%\n", percentage);

            //Cycles until the time is over
            while (inst->timelimit > seconds() - inst->start_time) {
                CPXsetdblparam(env, CPX_PARAM_TILIM, INTERNAL_TIME_LIMIT);
                double objval = INFINITY;
                int error =CPXgetobjval(env, lp, &objval);
                if (error) {
                    inst->best_value = 0;
                    return;
                    printf("FAILED TO RETRIEVE OBJVAL, error: %d\n", error);
                    print_error("Failed to retrieve the objval (CPXgetx() -> HARD_FIX)\n");
                }

                //Check if the solution just found is better than the previous one or it's value it's the same
                //as the previous cycle. If not I reduce the percentage of blocking the edges
                if (objval < inst->best_value) {
                    if (CPXgetx(env, lp, inst->solution, 0, inst->nvariables - 1)) {
                        print_error("Failed to optimize MIP (CPXgetx() -> HARD_FIX)\n");
                    }
                    inst->best_value = objval;
                } else {
                    if (percentage - 20 >= 0) {
                        percentage -= 20;
                    } else {
                        percentage = 0;
                    }
                    printf("->%d%%", percentage);
                }

                if (cnt != 0) {
                    //Restore the lower bounds to 0 only to the
                    //variables passed to the previously fixed to 1
                    for (int i = 0; i < cnt; i++) {
                        bd[i] = 0.0;
                    }

                    if (CPXchgbds(env, lp, cnt, indices, lu, bd)) {
                        print_error("Error on restoring variables - HARD_FIX");
                    }
                    //CPXwriteprob(env, lp, "variables_restored.lp", NULL);
                }


                //Iterate over each variable in the solution array
                cnt = 0;
                for (int i = 0; i < inst->nvariables; i++) {
                    if (inst->solution[i] > 0.5 && rand() % 100 < 0) {
                        indices[cnt] = i;
                        lu[cnt] = LOWER_BOUND;
                        bd[cnt++] = 1.0;
                    }
                }

                if (cnt != 0) {
                    //Change the bounds of the selected edges
                    if (CPXchgbds(env, lp, cnt, indices, lu, bd)) {
                        print_error("Error on restoring variables - HARD_FIX");
                    }
                }
                //CPXwriteprob(env, lp, "modified_variables.lp", NULL);

                CPXmipopt(env, lp);
            }
            printf("\nSOLUTION CALCULATED\n");
            printf("\nCOST: %f\n", inst->best_value);

            //Frees the memory used
            free(indices);
            free(lu);
            free(bd);
            break;

        case SOFT_FIX:
            //Code needed to say to cplex that we have a callback function
            if (CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE,
                                   sec_callback, inst)) {
                print_error("CPXcallbacksetfunc() error - SOFT_FIX");
            }


            inst->best_value = INFINITY;
            inst->start_time = seconds();
            //Change the internal time limit of cplex
            CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit / 3);

            //Variables used in the creation of the soft fixing constraints
            double k_opt = 2.0;
            double rhs = inst->nnodes - k_opt;
            char sense = GREAT_EQUAL;
            char **cname = (char **) calloc(1, sizeof(char *));
            cname[0] = (char *) calloc(100, sizeof(char));
            sprintf(cname[0], "SOFT_FIX_(K_OPT_%f)", k_opt);

            //Variables used to manage the new constraints
            int remove_row_flag = 0;
            int lastrow = CPXgetnumrows(env, lp);

            printf("First solving\n");
            CPXmipopt(env, lp);

            int cycle = 1;
            printf("Cycles->%d", cycle);
            //Cycles until the time is over
            while (inst->timelimit > seconds() - inst->start_time) {
                CPXsetdblparam(env, CPX_PARAM_TILIM, INTERNAL_TIME_LIMIT);
                double objval = INFINITY;
                int error =CPXgetobjval(env, lp, &objval);
                if (error) {
                    inst->best_value = 0;
                    return;
                    printf("FAILED TO RETRIEVE OBJVAL, error: %d\n", error);
                    print_error("Failed to retrieve the objval (CPXgetx() - SOFT_FIX)\n");
                }

                //Check if the solution just found is better than the previous one or it's value it's the same
                //as the previous cycle. If not change the k-opt value
                if (objval < inst->best_value) {
                    if (CPXgetx(env, lp, inst->solution, 0, inst->nvariables - 1)) {
                        print_error("Failed to optimize MIP (CPXgetx() - SOFT_FIX)\n");
                    }

                    inst->best_value = objval;

                } else {
                    k_opt += 2.0;
                    rhs = inst->nnodes - k_opt;
                    sprintf(cname[0], "SOFT_FIX_(K_OPT_%f)", k_opt);
                    printf("K increased, new value: %f\n", k_opt);
                }

                //If a constraint is been already added removes it
                if (remove_row_flag > 0) {
                    lastrow = CPXgetnumrows(env, lp);
                    if (CPXdelrows(env, lp, lastrow - 1, lastrow - 1)) {
                        print_error("Error while deleting the row (SOFT_FIX)");
                    }
                    //CPXwriteprob(env, lp, "logfiles/deleted_constraint(soft_fix).lp", NULL);
                }

                lastrow = CPXgetnumrows(env, lp);
                //Create a new constraint
                if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
                    print_error("Error while creating the constraints (SOFT_FIX)");
                }

                //Iterate over each variable in the solution array
                for (int i = 0; i < inst->nvariables; i++) {
                    if (inst->solution[i] > 0.9) {
                        if (CPXchgcoef(env, lp, lastrow, i, 1.0)) {
                            print_error("Error while changing coefs on SOFT_FIX");
                        }
                    }
                }

                //Flag for check if a constraint has been added
                remove_row_flag++;

                CPXmipopt(env, lp);
                printf("->%d", ++cycle);
            }
            printf("\n");
            //Delete the last constraint added in the last cycle
            lastrow = CPXgetnumrows(env, lp);
            if (CPXdelrows(env, lp, lastrow - 1, lastrow - 1)) {
                print_error("Error while deleting the last row (SOFT_FIX)");
            }

            printf("\nCOST: %f\n", inst->best_value);

            free(cname);
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

    //Create a new constraint for each component
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

// Callback switch
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

    free(solution_found);
    free(successors);
    free(component);

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

    if (ncomp == 1) {
        ccr_param inParam;
        inParam.context = context;
        inParam.inst = inst;
        if (CCcut_violated_cuts(inst->nnodes, ecount, elist, solution_found, 2.0 - EPS, create_cut_relaxation,
                                (void *) &inParam)) {
            print_error("Error on CCcut_violated_cuts()");
        }
    }

    free(compscount);
    free(comps);
    free(elist);
    free(solution_found);

    return 0;
}


int create_cut_relaxation(double cutval, int cutcount, int *cut, void *inParam) {
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
            rmatval[position++] = 1.0; //We can do it with memset
        }
    }

    if (CPXcallbackaddusercuts(parameters->context, 1, cutcount, &rhs, &sense, &rmatbeg, rmatind, rmatval,
                               &purgeable, &local)) {
        print_error("Error on CPXcallbackaddusercuts()");
    }

    free(rmatind);
    free(rmatval);

    return 0;
}


void greedy(instance *inst) {
    //Initialization of the array to save the best solution. The memory allocation is done in read_input.c
    for (int i = 0; i < inst->nnodes; i++) {
        inst->successors[i] = -1;
    }

    //Modify the way in which the solution is saved
    for (int starting_node = 0; starting_node < inst->nnodes; starting_node++) {
        //Initialization of the temporary array to save the solution
        int *tmp_successors = (int *) calloc(inst->nnodes, sizeof(int));
        for (int i = 0; i < inst->nnodes; i++) {
            tmp_successors[i] = -1;
        }

        int stop_flag = 0;

        int current_node = starting_node;
        double cost = 0;

        while (!stop_flag) {
            double min_distance = DBL_MAX;
            int candidate_node = -1;

            for (int i = 0; i < inst->nnodes; i++) {
                //If the node analyzed is the current ore the successors of i has been already defined
                //skips the iteration
                if (current_node == i || tmp_successors[i] != -1) continue;

                //If the node of the cycle is nearest in respect the previous one change this values
                if (distance(current_node, i, inst) < min_distance) {
                    min_distance = distance(current_node, i, inst);
                    candidate_node = i;
                }
            }

            if (candidate_node != -1) {
                tmp_successors[current_node] = candidate_node;
                current_node = candidate_node;
                cost += min_distance;
            } else {
                tmp_successors[current_node] = starting_node;
                cost += distance(current_node, starting_node, inst);
                stop_flag = 1;
            }
        }

        if (cost < inst->best_value) {
            inst->best_value = cost;
            for (int i = 0; i < inst->nnodes; i++) {
                inst->successors[i] = tmp_successors[i];
            }
        }

        free(tmp_successors);
    }
}


void extra_mileage(instance *inst) {
    //Initialization of the array to save the best solution. The memory allocation is done in read_input.c
    for (int i = 0; i < inst->nnodes; i++) {
        inst->successors[i] = -1;
    }

    //Find the farthest nodes in the problem
    int node1 = -1;
    int node2 = -1;
    double max_distance = 0;
    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = i + 1; j < inst->nnodes; j++) {
            if (i == j) continue;
            double tmp_distance = distance(i, j, inst);
            if (tmp_distance > max_distance) {
                max_distance = tmp_distance;
                node1 = i;
                node2 = j;
            }
        }
    }

    inst->successors[node1] = node2;
    inst->successors[node2] = node1;
    inst->best_value = 2 * max_distance;

    int start_node = node1;

    int node_added = 2;

    while (node_added < inst->nnodes) {

        int current_node = start_node;
        int candidate_node = -1;
        int candidate_start = -1;
        int candidate_end = -1;
        double min_distance = DBL_MAX;

        //Find the node which extra mileage is the smallest among the free nodes
        do {
            int next_node = inst->successors[current_node];

            for (int i = 0; i < inst->nnodes; i++) {
                if (i == current_node || i == next_node || inst->successors[i] != -1) continue;
                double delta = distance(current_node, i, inst) +
                               distance(i, next_node, inst) -
                               distance(current_node, next_node, inst);
                if (delta < min_distance) {
                    candidate_node = i;
                    candidate_start = current_node;
                    candidate_end = next_node;
                    min_distance = delta;
                }
            }

            current_node = next_node;
        } while (current_node != start_node);

        //Update the solution
        inst->successors[candidate_start] = candidate_node;
        inst->successors[candidate_node] = candidate_end;
        inst->best_value = inst->best_value -
                           distance(candidate_start, candidate_end, inst) +
                           distance(candidate_start, candidate_node, inst) +
                           distance(candidate_node, candidate_end, inst);

        node_added++;
    }
}


void two_opt_refining(instance *inst) {
    int flag_while = 0;

    while (!flag_while) {
        flag_while = 1;
        int candidate_i = -1;
        int candidate_j = -1;
        double max_delta = 0;
        for (int i = 0; i < inst->nnodes; i++) {
            for (int j = 0; j < inst->nnodes; j++) {
                //Skip rules
                if (inst->successors[i] == -1
                    || inst->successors[j] == -1
                    || i == j
                    || inst->successors[i] == j
                    || inst->successors[j] == i)
                    continue;

                //Difference of the 4 edges taken into considerations
                double delta = distance(i, inst->successors[i], inst) +
                               distance(j, inst->successors[j], inst) -
                               distance(i, j, inst) -
                               distance(inst->successors[i], inst->successors[j], inst);

                //If there is an improvement I save the values of the nodes to be changed
                if (delta > max_delta) {
                    candidate_i = i;
                    candidate_j = j;
                    max_delta = delta;
                    flag_while = 0;
                }
            }
        }

        //Change the order of successors
        if (max_delta > 0) {
            inst->best_value -= max_delta;
            int current = inst->successors[candidate_i];
            int previous = inst->successors[candidate_j];
            int end_node = previous;
            while (current != end_node) {
                int next = inst->successors[current];
                inst->successors[current] = previous;
                previous = current;
                current = next;
            }
            inst->successors[candidate_i] = candidate_j;
        }
    }
}


void three_opt_refining(instance *inst) {

    int flag_while = 0;

    while (!flag_while) {
        flag_while = 1;
        int flag_comb = 0;
        int candidate_i = -1;
        int candidate_j = -1;
        int candidate_k = -1;
        double max_delta = 0;

        for (int i = 0; i < inst->nnodes; i++) {
            int j = inst->successors[i];

            while (inst->successors[j] != i
                   && j != i) {

                //Skips the first assignment of the cycle
                if (j == inst->successors[i]) {
                    j = inst->successors[j];
                }

                int k = inst->successors[j];
                while (inst->successors[k] != i
                       && k != i) {

                    //Skips the first assignment of the cycle
                    if (k == inst->successors[j]) {
                        k = inst->successors[k];
                    }

                    //Previous value for the edges connected
                    double prev_value = distance(i, inst->successors[i], inst) +
                                        distance(j, inst->successors[j], inst) +
                                        distance(k, inst->successors[k], inst);

                    //All the combination of the new edges
                    double combs[4] = {
                            distance(i, j, inst) +
                            distance(inst->successors[i], k, inst) +
                            distance(inst->successors[j], inst->successors[k], inst),

                            distance(i, inst->successors[j], inst) +
                            distance(k, inst->successors[i], inst) +
                            distance(j, inst->successors[k], inst),

                            distance(i, inst->successors[j], inst) +
                            distance(k, j, inst) +
                            distance(inst->successors[i], inst->successors[k], inst),

                            distance(i, k, inst) +
                            distance(inst->successors[j], inst->successors[i], inst) +
                            distance(j, inst->successors[k], inst)
                    };

                    //Checks for new values of delta
                    for (int combo = 0; combo < 4; combo++) {
                        double delta = prev_value - combs[combo];
                        if (delta > max_delta) {
                            candidate_i = i;
                            candidate_j = j;
                            candidate_k = k;
                            max_delta = delta;
                            flag_comb = combo + 1;
                            flag_while = 0;
                        }
                    }

                    //Update k
                    k = inst->successors[k];
                }

                //Update j
                j = inst->successors[j];
            }

        }

        if (max_delta > 0) {
            inst->best_value -= max_delta;
            //Changes the order of the nodes following the best case
            switch (flag_comb) {
                case 1: {
                    inst->best_value -= max_delta;
                    int succ_j = inst->successors[candidate_j];
                    int succ_i = inst->successors[candidate_i];
                    int succ_k = inst->successors[candidate_k];

                    int current = succ_i;
                    int previous = candidate_k;
                    int end_node = succ_j;
                    while (current != end_node) {
                        int next = inst->successors[current];
                        inst->successors[current] = previous;
                        previous = current;
                        current = next;
                    }
                    inst->successors[candidate_i] = candidate_j;

                    current = succ_j;
                    previous = inst->successors[candidate_k];
                    end_node = succ_k;
                    while (current != end_node) {
                        int next = inst->successors[current];
                        inst->successors[current] = previous;
                        previous = current;
                        current = next;
                    }

                    break;
                }

                case 2: {
                    inst->best_value -= max_delta;
                    int succ_j = inst->successors[candidate_j];
                    int succ_i = inst->successors[candidate_i];
                    int succ_k = inst->successors[candidate_k];

                    inst->successors[candidate_i] = succ_j;
                    inst->successors[candidate_k] = succ_i;
                    inst->successors[candidate_j] = succ_k;
                    break;
                }

                case 3: {
                    inst->best_value -= max_delta;
                    int succ_j = inst->successors[candidate_j];
                    int succ_i = inst->successors[candidate_i];
                    int succ_k = inst->successors[candidate_k];


                    int current = succ_i;
                    int previous = succ_k;
                    int end_node = succ_j;
                    while (current != end_node) {
                        int next = inst->successors[current];
                        inst->successors[current] = previous;
                        previous = current;
                        current = next;
                    }
                    inst->successors[candidate_i] = succ_j;
                    inst->successors[candidate_k] = candidate_j;
                    break;
                }

                case 4: {
                    inst->best_value -= max_delta;
                    int succ_j = inst->successors[candidate_j];
                    int succ_i = inst->successors[candidate_i];
                    int succ_k = inst->successors[candidate_k];


                    int current = succ_j;
                    int previous = succ_i;
                    int end_node = succ_k;
                    while (current != end_node) {
                        int next = inst->successors[current];
                        inst->successors[current] = previous;
                        previous = current;
                        current = next;
                    }
                    inst->successors[candidate_j] = succ_k;
                    inst->successors[candidate_i] = candidate_k;
                    break;
                }

                default:
                    break;
            }
        }
    }
}


void vns(instance *inst) {
    double current_best = inst->best_value;
    int *best_successors = (int *) calloc(inst->nnodes, sizeof(int));
    for (int i = 0; i < inst->nnodes; i++) {
        best_successors[i] = inst->successors[i];
    }

    int flag = 0;
    while (!flag) {
        int cycles = 0;

        two_opt_refining(inst);
        // Updates the best solution
        if (inst->best_value < current_best) {
            current_best = inst->best_value;
            for (int i = 0; i < inst->nnodes; i++) {
                best_successors[i] = inst->successors[i];
            }
        } else {
            cycles++;
        }

        // Save the best solution in the instance and exit the node
        if (seconds() - inst->start_time > inst->timelimit) {
            inst->best_value = current_best;
            for (int i = 0; i < inst->nnodes; i++) {
                inst->successors[i] = best_successors[i];
            }
            flag = 1;
            continue;
        }

        // Performs kick
        if (cycles < 5) {
            three_kick_vns(inst);
        } else if (cycles >= 5 && cycles < 10) {
            five_kick_vns(inst);
        } else {
            seven_kick_vns(inst);
        }
    }
    printf("Final cost: %f\n", inst->best_value);
}


void three_kick_vns(instance *inst) {
    int first_node = -1;
    int second_node = -1;
    int third_node = -1;
    first_node = rand() % inst->nnodes;
    while (second_node == first_node || third_node == first_node ||
           second_node == third_node) {
        second_node = rand() % inst->nnodes;
        third_node = rand() % inst->nnodes;
    }

    // Order the elements
    int start_node = first_node;
    int curr = inst->successors[first_node];
    int nodes[3];
    int number = 1;
    nodes[0] = first_node;
    while (curr != start_node) {
        if (curr == second_node) {
            nodes[number++] = second_node;
        } else if (curr == third_node) {
            nodes[number++] = third_node;
        }
        if (number == 3) break;
        curr = inst->successors[curr];
    }

    int tmp_1 = inst->successors[nodes[0]];
    int tmp_2 = inst->successors[nodes[1]];
    int tmp_3 = inst->successors[nodes[2]];
    inst->successors[nodes[0]] = inst->successors[nodes[1]];
    inst->successors[nodes[1]] = tmp_3;
    inst->successors[nodes[2]] = tmp_1;

    inst->best_value = compute_solution_cost(inst, inst->successors);
}


void five_kick_vns(instance *inst) {
    int first_node = -1;
    int second_node = -1;
    int third_node = -1;
    int fourth_node = -1;
    int fifth_node = -1;
    first_node = rand() % inst->nnodes;
    while ((second_node = rand() % inst->nnodes) == first_node);
    while ((third_node = rand() % inst->nnodes) == first_node || third_node == second_node);
    while ((fourth_node = rand() % inst->nnodes) == first_node || fourth_node == second_node ||
           fourth_node == third_node);
    while ((fifth_node = rand() % inst->nnodes) == first_node || fifth_node == second_node ||
           fifth_node == third_node || fifth_node == fourth_node);

    // Order the elements
    int start_node = first_node;
    int curr = inst->successors[first_node];
    int nodes[5];
    int number = 1;
    nodes[0] = first_node;
    while (curr != start_node) {
        if (curr == second_node) {
            nodes[number++] = second_node;
        } else if (curr == third_node) {
            nodes[number++] = third_node;
        } else if (curr == fourth_node) {
            nodes[number++] = fourth_node;
        } else if (curr == fifth_node) {
            nodes[number++] = fifth_node;
        }
        if (number == 5) break;
        curr = inst->successors[curr];
    }

    int tmp_1 = inst->successors[nodes[0]];
    int tmp_2 = inst->successors[nodes[1]];
    int tmp_3 = inst->successors[nodes[2]];
    int tmp_4 = inst->successors[nodes[3]];
    int tmp_5 = inst->successors[nodes[4]];
    inst->successors[nodes[0]] = tmp_3;
    inst->successors[nodes[1]] = tmp_4;
    inst->successors[nodes[2]] = tmp_5;
    inst->successors[nodes[3]] = tmp_1;
    inst->successors[nodes[4]] = tmp_2;

    inst->best_value = compute_solution_cost(inst, inst->successors);
}


void seven_kick_vns(instance *inst) {
    int first_node = -1;
    int second_node = -1;
    int third_node = -1;
    int fourth_node = -1;
    int fifth_node = -1;
    int sixth_node = -1;
    int seventh_node = -1;
    first_node = rand() % inst->nnodes;
    while ((second_node = rand() % inst->nnodes) == first_node);
    while ((third_node = rand() % inst->nnodes) == first_node || third_node == second_node);
    while ((fourth_node = rand() % inst->nnodes) == first_node || fourth_node == second_node ||
           fourth_node == third_node);
    while ((fifth_node = rand() % inst->nnodes) == first_node || fifth_node == second_node ||
           fifth_node == third_node || fifth_node == fourth_node);
    while ((sixth_node = rand() % inst->nnodes) == first_node || sixth_node == second_node ||
           sixth_node == third_node || sixth_node == fourth_node || sixth_node == fifth_node);
    while ((seventh_node = rand() % inst->nnodes) == first_node || seventh_node == second_node ||
           seventh_node == third_node || seventh_node == fourth_node || seventh_node == fifth_node ||
           seventh_node == sixth_node);

    // Order the elements
    int start_node = first_node;
    int curr = inst->successors[first_node];
    int nodes[7];
    int number = 1;
    nodes[0] = first_node;
    while (curr != start_node) {
        if (curr == second_node) {
            nodes[number++] = second_node;
        } else if (curr == third_node) {
            nodes[number++] = third_node;
        } else if (curr == fourth_node) {
            nodes[number++] = fourth_node;
        } else if (curr == fifth_node) {
            nodes[number++] = fifth_node;
        } else if (curr == sixth_node) {
            nodes[number++] = sixth_node;
        } else if (curr == seventh_node) {
            nodes[number++] = seventh_node;
        }
        if (number == 7) break;
        curr = inst->successors[curr];
    }

    int tmp_1 = inst->successors[nodes[0]];
    int tmp_2 = inst->successors[nodes[1]];
    int tmp_3 = inst->successors[nodes[2]];
    int tmp_4 = inst->successors[nodes[3]];
    int tmp_5 = inst->successors[nodes[4]];
    int tmp_6 = inst->successors[nodes[5]];
    int tmp_7 = inst->successors[nodes[6]];
    inst->successors[nodes[0]] = tmp_4;
    inst->successors[nodes[1]] = tmp_5;
    inst->successors[nodes[2]] = tmp_6;
    inst->successors[nodes[3]] = tmp_1;
    inst->successors[nodes[4]] = tmp_2;
    inst->successors[nodes[5]] = tmp_7;
    inst->successors[nodes[6]] = tmp_3;

    inst->best_value = compute_solution_cost(inst, inst->successors);
}


void tabu_search(instance *inst) {
    double current_best = inst->best_value;
    int *best_successors = (int *) calloc(inst->nnodes, sizeof(int));
    for (int i = 0; i < inst->nnodes; i++) {
        best_successors[i] = inst->successors[i];
    }

    int *tabu_nodes = (int *) calloc(inst->nnodes, sizeof(int));
    for (int i = 0; i < inst->nnodes; i++) {
        tabu_nodes[i] = -1;
    }

    int iteration_counter = 0;
    // Should be dynamic (change over time)
    int tenure = inst->nnodes / 10;
    int flag = 0;
    while (!flag) {
        // Save the best solution in the instance and exit the node
        if (seconds() - inst->start_time > inst->timelimit) {
            inst->best_value = current_best;
            for (int i = 0; i < inst->nnodes; i++) {
                inst->successors[i] = best_successors[i];
            }
            flag = 1;
            continue;
        }


        int tabu_search = 0;
        // Performs tabu search
        while (!tabu_search) {
            int first_node;
            int second_node;
            double improvement;
            compute_bigger_cut(inst, &first_node, &second_node, &improvement);

            if ((tabu_nodes[first_node] != -1 && iteration_counter - tabu_nodes[first_node] <= tenure) ||
                (tabu_nodes[second_node] != -1 && iteration_counter - tabu_nodes[second_node] <= tenure)) {
                if (inst->best_value < current_best) {
                    current_best = inst->best_value;
                    for (int i = 0; i < inst->nnodes; i++) {
                        best_successors[i] = inst->successors[i];
                    }
                }
                tabu_search = 1;
                continue;
            }

            if (improvement == 0) {
                if (inst->best_value < current_best) {
                    current_best = inst->best_value;
                    for (int i = 0; i < inst->nnodes; i++) {
                        best_successors[i] = inst->successors[i];
                    }
                }
                tabu_search = 1;
                continue;
            }
            perform_cut(inst, first_node, second_node, improvement);
            iteration_counter++;
        }

        // 2-kick to exit the minima
        int first_node = rand() % inst->nnodes;
        int second_node;
        while ((second_node = rand() % inst->nnodes) == first_node || second_node == inst->successors[first_node]);
        double improvement = distance(first_node, inst->successors[first_node], inst) +
                             distance(second_node, inst->successors[second_node], inst) -
                             distance(first_node, second_node, inst) -
                             distance(inst->successors[first_node], inst->successors[second_node], inst);
        perform_cut(inst, first_node, second_node, improvement);

        //Adding the nodes to the tabu list
        tabu_nodes[first_node] = iteration_counter;
        tabu_nodes[second_node] = iteration_counter;

        iteration_counter++;
    }
    printf("Number of iterations: %d\n", iteration_counter);
    printf("Final cost: %f\n", inst->best_value);
}


void genetic(instance *inst, int population_size) {

    int **population = (int **) calloc(inst->nnodes, sizeof(int *));
    double *fitness = (double *) calloc(inst->nnodes, sizeof(double));

    double worst_fitness = -1;
    double average_fitness = -1;
    double initial_spread_fitness = -1;
    double best_fitness = INFINITY;
    int champion = -1;

    // Population generation
    for (int i = 0; i < population_size; i++) {
        population[i] = (int *) calloc(inst->nnodes, sizeof(int));
        fitness[i] = generate_random_solution(inst, population[i]);
        if (fitness[i] > worst_fitness) {
            worst_fitness = fitness[i];
        }
        if (fitness[i] < best_fitness) {
            best_fitness = fitness[i];
            champion = i;
        }
    }
    printf("Generated population\n");;

    int end_epochs = 0;
    int n_epochs = 0;
    while (!end_epochs) {

        //Generate new sons
        int number_of_children = population_size / 10;
        int **children = (int **) calloc(number_of_children, sizeof(int *));
        //Children generation
        for (int n_child = 0; n_child < number_of_children; n_child++) {
            int parent_1 = -1;
            int parent_2 = -1;

            while (parent_1 == -1) {
                int index = rand() % population_size;
                double normalized_fitness = (fitness[index] - best_fitness) / (worst_fitness - best_fitness);
                // Select from the top 50%
                if (normalized_fitness < 0.9) {
                    parent_1 = index;
                }
            }

            while (parent_2 == -1 || parent_2 == parent_1) {
                int index = rand() % population_size;
                double normalized_fitness = (fitness[index] - best_fitness) / (worst_fitness - best_fitness);
                if (normalized_fitness < 0.9) {
                    parent_2 = index;
                }
            }

            children[n_child] = (int *) calloc(inst->nnodes, sizeof(int));
            int cutting_point = inst->nnodes / 2;
            // First half chromosome
            for (int i = 0; i < cutting_point; i++) {
                children[n_child][i] = population[parent_1][i];
            }
            // Second half
            int next_node = cutting_point;
            int add_node = 1;
            for (int i = cutting_point; i < inst->nnodes; i++) {
                for (int j = 0; j < cutting_point; j++) {
                    if (population[parent_2][i] == children[n_child][j]) {
                        add_node = 0;
                        break;
                    }
                }
                if (add_node) {
                    children[n_child][next_node] = population[parent_2][i];
                    next_node++;
                } else {
                    add_node = 1;
                }
            }
            //Merging
            complete_merging(inst, children[n_child], next_node);
        }

        //Population killing and substituting with the children
        for (int n_child = 0; n_child < number_of_children; n_child++) {
            int killed = -1;

            while (killed == -1) {
                int index = rand() % population_size;
                double normalized_fitness = (fitness[index] - best_fitness) / (worst_fitness - best_fitness);
                // Select from the bottom 50%
                if (normalized_fitness > 0.2) {
                    killed = index;
                }
            }

            double initial1 = distance(children[n_child][0], children[n_child][inst->nnodes - 1], inst);
            for (int j = 0; j < inst->nnodes - 1; j++) {
                initial1 += distance(children[n_child][j], children[n_child][j + 1], inst);
            }

            fitness[killed] = distance(children[n_child][0], children[n_child][inst->nnodes - 1], inst);
            population[killed][inst->nnodes - 1] = children[n_child][inst->nnodes - 1];
            for (int i = 0; i < inst->nnodes - 1; i++) {
                population[killed][i] = children[n_child][i];
                fitness[killed] += distance(children[n_child][i], children[n_child][i + 1], inst);
            }

            free(children[n_child]);
        }
        free(children);

        // Updating fitness values
        double sum = 0.0;
        worst_fitness = -1;
        for (int i = 0; i < population_size; i++) {
            sum += fitness[i];
            if (fitness[i] < best_fitness) {
                best_fitness = fitness[i];
                champion = i;
            }
            if (fitness[i] > worst_fitness) {
                worst_fitness = fitness[i];
            }
        }
        average_fitness = sum / (double) population_size;
        double time = seconds();
        if (seconds() - inst->start_time > inst->timelimit) {
            inst->best_value = fitness[champion];
            list_to_successors(inst, population[champion], inst->successors, NULL);
            end_epochs = 1;
            continue;
        }

        if (n_epochs == 0) {
            initial_spread_fitness = worst_fitness - best_fitness;
        }

        // Perform mutations
        if ((worst_fitness - best_fitness) / initial_spread_fitness < 0.15) {
            int mutations = rand() % (inst->nnodes / 2);
            for (int i = 0; i < mutations; i++) {
                int index = -1;

                while (index == -1) {
                    int temp = rand() % population_size;
                    if (temp != champion) {
                        index = temp;
                    }
                }

                int first_node, second_node = -1;
                first_node = rand() % inst->nnodes;
                while ((second_node = rand() % inst->nnodes) == first_node);

                //Removing old edges from the fitness
                fitness[index] -= distance(population[index][first_node], population[index][first_node - 1], inst);
                fitness[index] -= distance(population[index][first_node], population[index][first_node + 1], inst);
                fitness[index] -= distance(population[index][second_node], population[index][second_node - 1], inst);
                fitness[index] -= distance(population[index][second_node], population[index][second_node + 1], inst);

                int temp = population[index][first_node];
                population[index][first_node] = population[index][second_node];
                population[index][second_node] = temp;

                //Adding the new edges to the fitness
                fitness[index] += distance(population[index][first_node], population[index][first_node - 1], inst);
                fitness[index] += distance(population[index][first_node], population[index][first_node + 1], inst);
                fitness[index] += distance(population[index][second_node], population[index][second_node - 1], inst);
                fitness[index] += distance(population[index][second_node], population[index][second_node + 1], inst);

                break;
            }
        }
        n_epochs++;
    }
}
