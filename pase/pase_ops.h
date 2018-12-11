/*
 * =====================================================================================
 *
 *       Filename:  gcge_para.h
 *
 *    Description:  
 *        基于GCGE的ops的单向量生成pase的单向量
 *        基于GCGE的ops的多向量生成pase的多向量
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

#ifndef  _PASE_OPS_H_
#define  _PASE_OPS_H_

#include "gcge_ops.h"
//把每一个操作都写好，这样以后进行算法设计的时候才能方便。
#if GCGE_USE_MPI
#include <mpi.h>
#endif


typedef struct PASE_OPS_ {

    GCGE_OPS *gcge_ops;

    void (*VecSetRandomValue)       (void *vec);
    void (*MatDotVec)               (void *Matrix, void *x, void *r);
    void (*VecAxpby)                (GCGE_DOUBLE a, void *x, GCGE_DOUBLE b, void *y); /* y = ax+by */
    void (*VecInnerProd)            (void *x, void *y, GCGE_DOUBLE *xTy);
    void (*VecLocalInnerProd)       (void *x, void *y, GCGE_DOUBLE *xTy);
    void (*VecCreateByVec)           (void *src_vec, void **des_vec);
    void (*VecCreateByMat)           (void *mat, void **vec);
    void (*VecDestroy)                 (void **vec);

    void (*MultiVecCreateByVec)      (void *vec, void ***multi_vec, GCGE_INT n_vec, 
                                     struct PASE_OPS_ *ops);
    void (*MultiVecCreateByMat)      (void *mat, void ***multi_vec, GCGE_INT n_vec, 
                                     struct PASE_OPS_ *ops);
    void (*MultiVecCreateByMultiVec) (void **init_vec, void ***multi_vec, GCGE_INT n_vec, 
                                     struct PASE_OPS_ *ops);
    void (*MultiVecDestroy)            (void ***MultiVec, GCGE_INT n_vec, struct PASE_OPS_ *ops);

    /* TODO */
    void (*MultiVecSetRandomValue)  (void **multi_vec, GCGE_INT start, GCGE_INT n_vec, struct PASE_OPS_ *ops);
    void (*MatDotMultiVec)          (void *mat, void **x, void **y, GCGE_INT *start, GCGE_INT *end, 
                                     struct PASE_OPS_ *ops);
    void (*MultiVecAxpby)           (GCGE_DOUBLE a, void **x, GCGE_DOUBLE b, void **y, 
                                     GCGE_INT *start, GCGE_INT *end, struct PASE_OPS_ *ops);
    void (*MultiVecAxpbyColumn)     (GCGE_DOUBLE a, void **x, GCGE_INT col_x, GCGE_DOUBLE b, 
                                     void **y, GCGE_INT col_y, struct PASE_OPS_ *ops);
    /* vec_y[j] = \sum_{i=sx}^{ex} vec_x[i] a[i-sx][j-sy] */
    void (*MultiVecLinearComb)      (void **x, void **y, GCGE_INT *start, GCGE_INT *end,
                                     GCGE_DOUBLE *a, GCGE_INT lda, 
                                     void *dmat, GCGE_INT lddmat, struct PASE_OPS_ *ops);
    void (*MultiVecInnerProd)       (void **V, void **W, GCGE_DOUBLE *a, char *is_sym, 
                                     GCGE_INT *start, GCGE_INT *end, GCGE_INT lda, struct PASE_OPS_ *ops);
    void (*MultiVecSwap)            (void **V_1, void **V_2, GCGE_INT *start, GCGE_INT *end, 
                                     struct PASE_OPS_ *ops);

    /* TODO kernal function should use this op to get j-th vector */
    void (*GetVecFromMultiVec)      (void **V, GCGE_INT j, void **x);
    void (*RestoreVecForMultiVec)   (void **V, GCGE_INT j, void **x);
    void (*PrintMultiVec)           (void **x, GCGE_INT n);
   
}PASE_OPS;
void PASE_OPS_Create(PASE_OPS **ops);
GCGE_INT PASE_OPS_Setup(PASE_OPS *ops);
void PASE_OPS_Free(PASE_OPS **ops);
#endif
