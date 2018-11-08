/*
 * =====================================================================================
 *
 *       Filename:  gcge_app_csr.h
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

#ifndef _GCGE_APP_CSR_H_
#define _GCGE_APP_CSR_H_

#include "gcge_type.h"
#include "gcge_ops.h"
#include "gcge_solver.h"

#include "csr.h"

void GCGE_CSR_SetOps(GCGE_OPS *ops);
void GCGE_SOLVER_SetCSROps(GCGE_SOLVER *solver);

GCGE_SOLVER *GCGE_CSR_Solver_Init(CSR_MAT *A, CSR_MAT *B, int num_eigenvalues, int argc, char* argv[]);
#endif
