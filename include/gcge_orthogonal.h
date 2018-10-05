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


#include "gcge_type.h"

#include "gcge_para.h"
#include "gcge_ops.h"

#include "gcge_workspace.h"

//GCGE算法中用到的正交化操作
void GCGE_ClassicalOrthogonalize(void **V, GCGE_INT start, GCGE_INT *end, 
                                 void *B, GCGE_OPS *ops, GCGE_ORTH_PARA *orth_para, 
                                 void **V_tmp, GCGE_DOUBLE *d_tmp);
void GCGE_ModifiedOrthogonalize(void **V, GCGE_INT start, GCGE_INT *end, 
                                void *B, GCGE_OPS *ops, GCGE_ORTH_PARA *orth_para, 
                                void **V_tmp, GCGE_DOUBLE *d_tmp);
void GCGE_BlockedOrthogonalize(void **V, GCGE_INT start, GCGE_INT *end, 
                               void *B, GCGE_OPS *ops, GCGE_ORTH_PARA *orth_para, 
                               void **V_tmp, GCGE_DOUBLE *d_tmp);

void GCGE_ModifiedOrthogonalInSubspace(double *V, GCGE_INT ldV, GCGE_INT start, GCGE_INT *end, 
                                       void *B, GCGE_INT ldB, GCGE_ORTH_PARA *orth_para);
#endif
