/*
 * =====================================================================================
 *
 *       Filename:  gcge_solver.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/27/2018 10:03:34 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef _GCGE_SOLVER_H_
#define _GCGE_SOLVER_H_

#include "gcge_type.h"
#include "gcge_ops.h"
#include "gcge_para.h"
#include "gcge_workspace.h"
#include "gcge_eigsol.h"

typedef struct GCGE_SOLVER_ {

   void *A;
   void *B;
   void **evec;
   GCGE_DOUBLE *eval;

   GCGE_PARA *para;
   GCGE_OPS  *ops;
   GCGE_WORKSPACE *workspace;
   
}GCGE_SOLVER;


GCGE_INT GCGE_SOLVER_Create(GCGE_SOLVER **solver, GCGE_INT argc, char* argv[]);
void GCGE_SOLVER_Free(GCGE_SOLVER **solver);
void GCGE_SOLVER_Free_All(GCGE_SOLVER **solver);

void GCGE_SOLVER_SetMatA(GCGE_SOLVER *solver, void *A);
void GCGE_SOLVER_SetMatB(GCGE_SOLVER *solver, void *B);
void GCGE_SOLVER_SetEigenvalues(GCGE_SOLVER *solver, GCGE_DOUBLE *eval);
void GCGE_SOLVER_SetEigenvectors(GCGE_SOLVER *solver, void **evec);
void GCGE_SOLVER_SetNumEigen(GCGE_SOLVER *solver, GCGE_INT nev);
void GCGE_SOLVER_Setup(GCGE_SOLVER *solver);
void GCGE_SOLVER_SetOpsLinearSolverWorkspace(GCGE_SOLVER *solver, void *linear_solver_workspace);
//从solver中取出特征对
void GCGE_SOLVER_GetEigenvalues(GCGE_SOLVER *solver, GCGE_DOUBLE **eval);
void GCGE_SOLVER_GetEigenvectors(GCGE_SOLVER *solver, void ***evec);
//求解特征值的子程序
void GCGE_SOLVER_Solve(GCGE_SOLVER *solver);

#endif
