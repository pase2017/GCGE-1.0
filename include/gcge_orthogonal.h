/*
 * =====================================================================================
 *
 *       Filename:  gcge_orthogonal.h
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
#ifndef _GCGE_ORTHOGONAL_H_
#define _GCGE_ORTHOGONAL_H_

#include <string.h>
#include "gcge_type.h"

#include "gcge_para.h"
#include "gcge_ops.h"

#include "gcge_workspace.h"
#include "gcge_utilities.h"
#include "gcge_matvec.h"

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
void GCGE_Orthogonal(void **V, GCGE_INT start, GCGE_INT *end, 
      void *B, GCGE_OPS *ops, GCGE_PARA *para, 
      void **V_tmp, GCGE_DOUBLE *d_tmp);

void GCGE_OrthogonalSubspace(double *V, GCGE_INT ldV, GCGE_INT nrows, GCGE_INT start, GCGE_INT *end, 
      void *B, GCGE_INT ldB, GCGE_ORTH_PARA *orth_para);
void GCGE_BOrthogonal(void **V, GCGE_INT start, GCGE_INT *end, 
                      void *B, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace);
void GCGE_CBOrthogonal(void **V, GCGE_INT start, GCGE_INT *end, 
                      void *B, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace);
void GCGE_SCBOrthogonal(void **V, GCGE_INT start, GCGE_INT *end, 
                      void *B, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace);
#endif
