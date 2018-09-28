/*
 * =====================================================================================
 *
 *       Filename:  gcge_solver.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2018年09月24日 09时50分16秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */


#include <stdio.h>
#include <stdlib.h>
#include "gcge_solver.h"


GCGE_INT GCGE_SOLVER_Create(GCGE_SOLVER **solver, GCGE_INT argc, char* argv[])
{
    *solver = (GCGE_SOLVER*)malloc(sizeof(GCGE_SOLVER));
    GCGE_PARA_Create(&(*solver)->para);
    GCGE_OPS_Create(&(*solver)->ops);
    GCGE_INT error = GCGE_PARA_SetFromCommandLine((*solver)->para, argc, argv);
    /* TODO should set NULL and Setup to modify */
    GCGE_WORKSPACE_Create(&(*solver)->workspace);
    (*solver)->A = NULL;
    (*solver)->B = NULL;
    (*solver)->eval = NULL;
    (*solver)->evec = NULL;
    return error;
}
void GCGE_SOLVER_Free  (GCGE_SOLVER **solver)
{
    GCGE_WORKSPACE_Free(&(*solver)->workspace, (*solver)->para, (*solver)->ops);
    GCGE_PARA_Free(&(*solver)->para);
    GCGE_OPS_Free(&(*solver)->ops);
    free(*solver); *solver = NULL;
}
void GCGE_SOLVER_Setup(GCGE_SOLVER *solver)
{
    GCGE_PARA_Setup(solver->para);
    GCGE_OPS_Setup (solver->ops);
    GCGE_WORKSPACE_Setup(solver->workspace, solver->para, solver->ops, solver->A);
}

void GCGE_SOLVER_SetMatA(GCGE_SOLVER *solver, void *A)
{
    solver->A = A;
}
void GCGE_SOLVER_SetMatB(GCGE_SOLVER *solver, void *B)
{
    solver->B = B;
}
void GCGE_SOLVER_SetEigenvalues(GCGE_SOLVER *solver, GCGE_DOUBLE *eval)
{
    solver->eval = eval;
}
void GCGE_SOLVER_SetEigenvectors(GCGE_SOLVER *solver, void **evec)
{
    solver->evec = evec;
}
void GCGE_SOLVER_SetNumEigen(GCGE_SOLVER *solver, GCGE_INT nev)
{
   GCGE_PARA_SetNumEigen(solver->para, nev);
}

void GCGE_SOLVER_Solve(GCGE_SOLVER *solver)
{
    GCGE_EigenSolver(solver->A, solver->B, solver->eval, solver->evec, solver->para, solver->ops, solver->workspace);
}
