/*
 * =====================================================================================
 *
 *       Filename:  gcge_app_hypre.h
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

#ifndef _GCGE_APP_HYPRE_H_
#define _GCGE_APP_HYPRE_H_

#include "gcge_config.h"
#include "gcge_ops.h"
#include "gcge_solver.h"

#include "_hypre_parcsr_mv.h"
#include "_hypre_utilities.h"

#include "HYPRE_parcsr_ls.h"


void GCGE_HYPRE_SetOps(GCGE_OPS *ops);
void GCGE_SOLVER_SetHYPREOps(GCGE_SOLVER *solver);

GCGE_SOLVER *GCGE_HYPRE_Solver_Init(HYPRE_ParCSRMatrix A, HYPRE_ParCSRMatrix B, int num_eigenvalues, int argc, char* argv[]);
GCGE_SOLVER* GCGE_HYPRE_Solver_Init_LinearSolverGivenByUser(HYPRE_ParCSRMatrix A, HYPRE_ParCSRMatrix B, 
      HYPRE_Solver linear_solver, int ls_id, int num_eigenvalues, int argc, char* argv[]);
#endif
