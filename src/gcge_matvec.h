/*
 * =====================================================================================
 *
 *       Filename:  gcge_matvec.h
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

#ifndef _GCGE_MATVEC_H_
#define _GCGE_MATVEC_H_
#include "gcge_config.h"
#include "gcge_ops.h"

//GCGE算法中用到的一些线性代数操作
GCGE_DOUBLE GCGE_VecMatrixVec(void *a, void *Matrix, void *b, void *temp, GCGE_OPS *ops);
GCGE_DOUBLE GCGE_VecNorm(void *x, GCGE_OPS *ops);

//子空间线性代数操作
GCGE_DOUBLE GCGE_ArrayDotArrayInSubspace(GCGE_DOUBLE *a, GCGE_DOUBLE *b, GCGE_INT length);
GCGE_DOUBLE GCGE_ArrayNormInSubspace(GCGE_DOUBLE *a, GCGE_INT length);
void GCGE_ArrayAXPBYInSubspace(GCGE_DOUBLE a, GCGE_DOUBLE *x, GCGE_DOUBLE b, GCGE_DOUBLE *y, 
        GCGE_INT length);
void GCGE_ArrayCopyInSubspace(GCGE_DOUBLE *x, GCGE_DOUBLE *y, GCGE_INT length);
void GCGE_ArrayScaleInSubspace(GCGE_DOUBLE alpha, GCGE_DOUBLE *a, GCGE_INT length);
#endif
