/*
 * =====================================================================================
 *
 *       Filename:  gcge_app_phg.h
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

#ifndef _GCGE_APP_PHG_H_
#define _GCGE_APP_PHG_H_

#include "gcge_type.h"
#include "gcge_ops.h"
#include "gcge_solver.h"

#include "phg.h"


void GCGE_PHG_SetOps(GCGE_OPS *ops);
void GCGE_SOLVER_SetPHGOps(GCGE_SOLVER *solver);

GCGE_SOLVER *GCGE_PHG_Solver_Init(MAT* A, MAT* B, int num_eigenvalues, int argc, char* argv[]);
GCGE_SOLVER* GCGE_PHG_Solver_Init_LinearSolverGivenByUser(MAT* A, MAT* B, 
      SOLVER* linear_solver, int num_eigenvalues, int argc, char* argv[]);
#endif
