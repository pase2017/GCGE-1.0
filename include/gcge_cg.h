/*
 * =====================================================================================
 *
 *       Filename:  gcge_cg.h
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
#ifndef _GCGE_CG_H_
#define _GCGE_CG_H_

#include "gcge_type.h"

#include "gcge_para.h"
#include "gcge_ops.h"

#include "gcge_workspace.h"
#include "gcge_matvec.h"

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
void GCGE_CG(void *Matrix, void *b, void *x, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace);
#endif
