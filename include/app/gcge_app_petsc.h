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

//#include <petscksp.h>
//#include <petscoptions.h>
//#include <petsc/private/vecimpl.h>
//#include <petscblaslapack.h>
#include <petscvec.h>
#include <petscmat.h>

void GCGE_PETSC_SetOps(GCGE_OPS *ops);
void GCGE_SOLVER_SetPETSCOps(GCGE_SOLVER *solver);

GCGE_SOLVER *GCGE_PETSC_Solver_Init(Mat A, Mat B, int num_eigenvalues, int argc, char* argv[]);
#endif
