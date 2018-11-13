/*
 * =====================================================================================
 *
 *       Filename:  gcge_eigsol.h
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
#ifndef _GCGE_EIGSOL_H_
#define _GCGE_EIGSOL_H_


#include "gcge_type.h"

#include "gcge_utilities.h"
#include "gcge_statistics.h"

#include "gcge_para.h"
#include "gcge_ops.h"

#include "gcge_workspace.h"
#include "gcge_matvec.h"

#include "gcge_orthogonal.h"
#include "gcge_rayleighritz.h"
#include "gcge_xpw.h"

//GCGE算法中用到的正交化操作
/**
 * @brief 
 *
 * @param V
 * @param start
 * @param[in|out] end
 * @param B
 * @param ops
 * @param orth_para
 * @param workspace
 */
void GCGE_EigenSolver(void *A, void *B, GCGE_DOUBLE *eval, void **evec, 
        GCGE_PARA *para, GCGE_OPS *ops, GCGE_WORKSPACE *workspace);
void GCGE_CheckConvergence(void *A, void *B, GCGE_DOUBLE *eval, void **evec, 
        GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace);
void GCGE_CheckEvalMultiplicity(GCGE_INT start, GCGE_INT nev, GCGE_INT end, GCGE_INT *dim_x, GCGE_DOUBLE *eval);

#endif
