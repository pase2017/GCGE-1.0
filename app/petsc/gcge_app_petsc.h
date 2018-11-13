/*
 * =====================================================================================
 *
 *       Filename:  gcge_app_petsc.h
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

#ifndef _GCGE_APP_PETSC_H_
#define _GCGE_APP_PETSC_H_

#include "gcge_type.h"
#include "gcge_ops.h"
#include "gcge_solver.h"

//#include <petscoptions.h>
//#include <petsc/private/vecimpl.h>
//#include <petscblaslapack.h>
#include <petscoptions.h>
#include <petscviewer.h>
#include <petscsys.h>

#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#include <petscpc.h>

void PETSC_ReadMatrixBinary(Mat *A, const char *filename);
void PETSC_LinearSolverCreate(KSP *ksp, Mat A, Mat T);
void PETSC_VecLocalInnerProd(Vec x, Vec y, double *value);

void GCGE_PETSC_SetOps(GCGE_OPS *ops);
void GCGE_SOLVER_SetPETSCOps(GCGE_SOLVER *solver);

GCGE_SOLVER *GCGE_PETSC_Solver_Init(Mat A, Mat B, int num_eigenvalues, int argc, char* argv[]);
GCGE_SOLVER* GCGE_PETSC_Solver_Init_KSPDefault(Mat A, Mat B, Mat P, int num_eigenvalues, int argc, char* argv[]);
GCGE_SOLVER* GCGE_PETSC_Solver_Init_KSPGivenByUser(Mat A, Mat B, KSP ksp, int num_eigenvalues, int argc, char* argv[]);
//设置ksp为线性解法器
void GCGE_SOLVER_SetPETSCOpsLinearSolver(GCGE_SOLVER *solver, KSP ksp);
#endif
