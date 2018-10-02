/*
 * =====================================================================================
 *
 *       Filename:  gcge_matvec.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2018年09月24日 09时57分13秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gcge_matvec.h"

//计算U = V * x
//由一组基向量组V与子空间向量x相乘得到一个长向量U
//V是用来计算基底向量组(GCGE_Vec),U是计算得到的向量(GCGE_Vec)
//void GCGE_SumSeveralVecs(void **V, GCGE_DOUBLE *x, void *U, GCGE_INT n_vec, 
//void GCGE_VecsLinearCombination(void **V, GCGE_DOUBLE *x, void *U, GCGE_INT n_vec, GCGE_OPS *ops)
//{
//    GCGE_INT i;
//    //U = 0
//    ops->VecAxpby(0.0, U, 0.0, U);
//    for( i=0; i<n_vec; i++ )
//    {
//        //U += x[i] * V[i]
//        ops->VecAxpby(x[i], V[i], 1.0, U);
//    }
//}

//计算a^T*Matrix*b,两个向量的内积
//a,b表示要计算内积的两个向量(GCGE_Vec)，
//Matrix表示计算内积用到的矩阵(GCGE_MATRIX),
//temp是临时向量存储空间(GCGE_Vec)
//GCGE_DOUBLE GCGE_VecMatrixVec(void *a, void *Matrix, void *b, void *temp, GCGE_OPS *ops)
//{
//    GCGE_DOUBLE value=0.0;
//    ops->MatDotVec(Matrix, b, temp);
//    ops->VecInnerProd(a, temp, &value);
//    return value;
//}

//计算向量范数
//x是向量(GCGE_Vec)
//GCGE_DOUBLE GCGE_VecNorm(void *x, GCGE_OPS *ops)
//{
//    GCGE_DOUBLE xTx;
//    ops->VecInnerProd(x, x, &xTx);
//    return sqrt(xTx);
//}

//子空间向量操作
//计算向量a和b的内积
GCGE_DOUBLE GCGE_VecDotVecSubspace(GCGE_DOUBLE *a, GCGE_DOUBLE *b, GCGE_INT n)
{
#if 1
    GCGE_INT i;
    GCGE_DOUBLE value = 0.0;
    for( i=0; i<n; i++ )
        value += a[i]*b[i];
    return value;
#else
    GCGE_INT inc = 1;
    return BLASname(dot)(&n, a, &inc, b, &inc);
#endif
}

//计算向量a的范数
GCGE_DOUBLE GCGE_VecNormSubspace(GCGE_DOUBLE *a, GCGE_INT n)
{
#if 1
    GCGE_INT i;
    GCGE_DOUBLE value = 0.0;
    for( i=0; i<n; i++ )
        value += a[i]*a[i];
    return sqrt(value);
#else
    GCGE_INT inc = 1;
    return BLASname(nrm2)(&n, a, &inc);
#endif
}

//计算y=a*x+b*y
void GCGE_VecAXPBYSubspace(GCGE_DOUBLE a, GCGE_DOUBLE *x, GCGE_DOUBLE b, GCGE_DOUBLE *y, 
        GCGE_INT n)
{
#if 1
    GCGE_INT i;
    for( i=0; i<n; i++ )
        y[i] = a*x[i] + b*y[i];
#else
    GCGE_INT inc = 1;
    if(b != 1.0)
    {
        BLASname(scal)(&n, &b, y, &inc);
    }
    BLASname(axpy)(&n, &a, x, &inc, y, &inc);
#endif
}

/* TODO y = x, should shift parameters */
void GCGE_VecCopySubspace(GCGE_DOUBLE *x, GCGE_DOUBLE *y, GCGE_INT n)
{
#if 1
    memcpy(y, x, n*sizeof(GCGE_DOUBLE));
#else
    GCGE_INT inc = 1;
    BLASname(copy)(&n, x, &inc, y, &inc);
#endif
}

//a=alpha*a
void GCGE_VecScaleSubspace(GCGE_DOUBLE alpha, GCGE_DOUBLE *a, GCGE_INT n)
{
#if 1
    GCGE_INT i;
    for( i=0; i<n; i++ )
        a[i] *= alpha;
#else
    GCGE_INT inc = 1;
    BLASname(scal)(&n, &alpha, a, &inc);
#endif
}

