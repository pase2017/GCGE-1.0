/*
 * =====================================================================================
 *
 *       Filename:  gcge_xpw.h
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
#ifndef _GCGE_XPW_H_
#define _GCGE_XPW_H_


#include "gcge_type.h"
#include "gcge_utilities.h"
#include "gcge_statistics.h"

#include "gcge_para.h"
#include "gcge_ops.h"

#include "gcge_workspace.h"
#include "gcge_matvec.h"

#include "gcge_cg.h"
#include "gcge_orthogonal.h"

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

void GCGE_ComputeX(void **V, GCGE_DOUBLE *subspace_evec, void **X, 
        GCGE_INT size_X, GCGE_INT size_V,
        GCGE_OPS *ops, GCGE_WORKSPACE *workspace);

void GCGE_ComputeW(void *A, void *B, void **V, GCGE_DOUBLE *eval, 
        GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace);

void GCGE_ComputeP(GCGE_DOUBLE *subspace_evec, void **V, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace);

#endif
