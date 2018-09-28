/*
 * =====================================================================================
 *
 *       Filename:  gcge_para.h
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

#ifndef  _GCGE_OPS_H_
#define  _GCGE_OPS_H_

#include "gcge_type.h"


typedef struct GCGE_OPS_ {

    void (*VecSetRandomValue)       (void *vec);
    void (*MatDotVec)               (void *Matrix, void *x, void *r);
    void (*VecAxpby)                (GCGE_DOUBLE a, void *x, GCGE_DOUBLE b, void *y); /* y = ax+by */
    void (*VecInnerProd)            (void *x, void *y, GCGE_DOUBLE *xTy);
    void (*VecLocalInnerProd)       (void *x, void *y, GCGE_DOUBLE *xTy);
    void (*BuildVecByVec)           (void *src_vec, void **des_vec);
    void (*BuildVecByMat)           (void *mat, void **vec);
    void (*FreeVec)                 (void **vec);


    /* option */
    /* TODO */
    void (*LinearSolver)            (void *Matrix, void *b, void *x, struct GCGE_OPS_ *ops);
    void *linear_solver_workspace;

    /* DenseMatCreate, DenseMatDestroy should in function Orthogonal 
     * Add struct member name void *orth_workspace to save tmp variables */
    void (*Orthogonal)              (void **V, GCGE_INT start, GCGE_INT *end, 
                                     void *B, GCGE_DOUBLE orth_zero_tol, 
                                     void **work); /* TODO */
    void (*DenseMatCreate)          (void **densemat, GCGE_INT nrows, GCGE_INT ncols);
    void (*DenseMatDestroy)         (void **mat);
 

    /* DEEP */

    void (*BuildMultiVecByVec)      (void *vec, void ***multi_vec, GCGE_INT n_vec, 
                                     struct GCGE_OPS_ *ops);
    void (*BuildMultiVecByMat)      (void *mat, void ***multi_vec, GCGE_INT n_vec, 
                                     struct GCGE_OPS_ *ops);
    void (*BuildMultiVecByMultiVec) (void **init_vec, void ***multi_vec, GCGE_INT n_vec, 
                                     struct GCGE_OPS_ *ops);
    void (*FreeMultiVec)            (void ***MultiVec, GCGE_INT n_vec, struct GCGE_OPS_ *ops);

    /* TODO */
    void (*MultiVecSetRandomValue)  (void **multi_vec, GCGE_INT n_vec, struct GCGE_OPS_ *ops);
    void (*MatDotMultiVec)          (void *mat, void **x, void **y, GCGE_INT *start, GCGE_INT *end, 
                                     struct GCGE_OPS_ *ops);
    void (*MultiVecAxpby)           (GCGE_DOUBLE a, void **x, GCGE_DOUBLE b, void **y, 
                                     GCGE_INT *start, GCGE_INT *end, struct GCGE_OPS_ *ops);
    void (*MultiVecAxpbyColumn)     (GCGE_DOUBLE a, void **x, GCGE_INT col_x, GCGE_DOUBLE b, 
                                     void **y, GCGE_INT col_y, struct GCGE_OPS_ *ops);
    /* vec_y[j] = \sum_{i=sx}^{ex} vec_x[i] a[i][j] */
    void (*MultiVecLinearComb)      (void **x, void **y, GCGE_INT *start, GCGE_INT *end,
                                     GCGE_DOUBLE *a, GCGE_INT lda, 
                                     void *dmat, GCGE_INT lddmat, struct GCGE_OPS_ *ops);
    void (*MultiVecInnerProd)       (void **V, void **W, GCGE_DOUBLE *a, char *is_sym, 
                                     GCGE_INT *start, GCGE_INT *end, GCGE_INT lda, struct GCGE_OPS_ *ops);
    void (*MultiVecSwap)            (void **V_1, void **V_2, GCGE_INT *start, GCGE_INT *end, 
                                     struct GCGE_OPS_ *ops);

    /* TODO kernal function should use this op to get j-th vector */
    void (*GetVecFromMultiVec)      (void **V, GCGE_INT j, void **x);
    void (*RestoreVecForMultiVec)   (void **V, GCGE_INT j, void **x);


    /* DEEP option */
    //DenseMatEigenSolver可取lapack中的dsyev,dsyevx,dsyevr,dsyevd
    //对dsyevx:
    /* param:
     * jobz: "V"，表示计算特征值及特征向量
     * range: "I"，表示计算第il到iu个特征对
     * uplo: "U"，表示给的是对称矩阵对上三角部分
     * nrows: 矩阵A的行数
     * a: 矩阵A
     * lda: 矩阵A的leading dimension
     * vl: 不用
     * vu: 不用
     * il: 个数的起始位置
     * iu: 个数的终点位置
     * abstol:
     * nev: 特征值个数
     * eval: 特征值
     * evec: 特征向量
     * lde: evec的leading dimension
     * isuppz: -----不要设，跳过这个参数------
     * work, 
     * lwork: 8*N
     * iwork, 
     * liwork: -----不要设，跳过这个参数------
     * ifail: int*, 长度2*nev
     * info
     */
    //对dsyevr:
    /* param:
     * jobz: "V"，表示计算特征值及特征向量
     * range: "I"，表示计算第il到iu个特征对
     * uplo: "U"，表示给的是对称矩阵对上三角部分
     * nrows: 矩阵A的行数
     * a: 矩阵A
     * lda: 矩阵A的leading dimension
     * vl: 不用
     * vu: 不用
     * il: 个数的起始位置
     * iu: 个数的终点位置
     * abstol:
     * nev: 特征值个数
     * eval: 特征值
     * evec: 特征向量
     * lde: evec的leading dimension
     * isuppz: int *
     * work, 
     * lwork: 26*N
     * iwork, 
     * liwork: 10*N
     * ifail: -----不要设，跳过这个参数------
     * info
     */
    //对dsyev:
    /* param:
     * jobz: "V"，表示计算特征值及特征向量
     * range: -----不要设，跳过这个参数------
     * uplo: "U"，表示给的是对称矩阵对上三角部分
     * nrows: 矩阵A的行数
     * a: 矩阵A
     * lda: 矩阵A的leading dimension
     * vl: -----不要设，跳过这个参数------
     * vu: -----不要设，跳过这个参数------
     * il: -----不要设，跳过这个参数------
     * iu: -----不要设，跳过这个参数------
     * abstol:-----不要设，跳过这个参数------
     * nev: -----不要设，跳过这个参数------
     * eval: 特征值
     * evec: -----不要设，跳过这个参数------
     * lde: -----不要设，跳过这个参数------
     * isuppz: -----不要设，跳过这个参数------
     * work, 
     * lwork: 3*N
     * iwork, -----不要设，跳过这个参数------
     * liwork: -----不要设，跳过这个参数------
     * ifail: -----不要设，跳过这个参数------
     * info
     */
    //对dsyevd:
    /* param:
     * jobz: "V"，表示计算特征值及特征向量
     * range: -----不要设，跳过这个参数------
     * uplo: "U"，表示给的是对称矩阵对上三角部分
     * nrows: 矩阵A的行数
     * a: 矩阵A
     * lda: 矩阵A的leading dimension
     * vl: -----不要设，跳过这个参数------
     * vu: -----不要设，跳过这个参数------
     * il: -----不要设，跳过这个参数------
     * iu: -----不要设，跳过这个参数------
     * abstol:-----不要设，跳过这个参数------
     * nev: -----不要设，跳过这个参数------
     * eval: 特征值
     * evec: -----不要设，跳过这个参数------
     * lde: -----不要设，跳过这个参数------
     * isuppz: -----不要设，跳过这个参数------
     * work, 
     * lwork: 1+6*N+2*N**2
     * iwork, 
     * liwork: 3+5*N
     * ifail: -----不要设，跳过这个参数------
     * info
     */
    void (*DenseMatEigenSolver)     (char *jobz, char *range, char *uplo, 
                                     GCGE_INT *nrows, GCGE_DOUBLE *a, GCGE_INT *lda, 
                                     GCGE_DOUBLE *vl, GCGE_DOUBLE *vu, GCGE_INT *il, GCGE_INT *iu, 
                                     GCGE_DOUBLE *abstol, GCGE_INT *nev, 
                                     GCGE_DOUBLE *eval, GCGE_DOUBLE *evec, GCGE_INT *lde, 
                                     GCGE_INT *isuppz, GCGE_DOUBLE *work, GCGE_INT *lwork, 
                                     GCGE_INT *iwork, GCGE_INT *liwork, 
                                     GCGE_INT *ifail, GCGE_INT *info);
    void (*DenseMatDotDenseMat)     (char *transa, char *transb, GCGE_INT *nrows, GCGE_INT *ncols, 
                                     GCGE_INT *mid, GCGE_DOUBLE *alpha, GCGE_DOUBLE *a,
                                     GCGE_INT *lda, GCGE_DOUBLE *b, GCGE_INT *ldb, 
                                     GCGE_DOUBLE *beta, GCGE_DOUBLE *c, GCGE_INT *ldc);
    void (*DenseSymMatDotDenseMat)  (char *side, char *uplo, GCGE_INT *nrows, GCGE_INT *ncols,
                                     GCGE_DOUBLE *alpha, GCGE_DOUBLE *a, GCGE_INT *lda,
                                     GCGE_DOUBLE *b, GCGE_INT *ldb, GCGE_DOUBLE *beta,
                                     GCGE_DOUBLE *c, GCGE_INT *ldc);

    void (*PrintMultiVec)           (void **x, GCGE_INT n);
   
}GCGE_OPS;

extern void dsyevx_(char *jobz, char *range, char *uplo, 
        GCGE_INT    *nrows,  GCGE_DOUBLE *a,    GCGE_INT *lda, 
        GCGE_DOUBLE *vl,     GCGE_DOUBLE *vu,   GCGE_INT *il,  GCGE_INT *iu, 
        GCGE_DOUBLE *abstol, GCGE_INT *nev, 
        GCGE_DOUBLE *eval,   GCGE_DOUBLE *evec, GCGE_INT *lde, 
        GCGE_DOUBLE *work,   GCGE_INT *lwork, 
        GCGE_INT    *iwork,  GCGE_INT *liwork,  GCGE_INT *info);

extern void dgemm_(char *transa, char *transb, 
        GCGE_INT    *nrows, GCGE_INT    *ncols, GCGE_INT    *mid, 
        GCGE_DOUBLE *alpha, GCGE_DOUBLE *a,     GCGE_INT    *lda, 
        GCGE_DOUBLE *b,     GCGE_INT    *ldb,   GCGE_DOUBLE *beta, 
        GCGE_DOUBLE *c,     GCGE_INT    *ldc);

extern void dsymm_(char *side, char *uplo, 
        GCGE_INT    *nrows, GCGE_INT    *ncols,
        GCGE_DOUBLE *alpha, GCGE_DOUBLE *a,   GCGE_INT    *lda,
        GCGE_DOUBLE *b,     GCGE_INT    *ldb, GCGE_DOUBLE *beta,
        GCGE_DOUBLE *c,     GCGE_INT    *ldc);


void GCGE_OPS_Create(GCGE_OPS **ops);
GCGE_INT GCGE_OPS_Setup(GCGE_OPS *ops);
void GCGE_OPS_Free(GCGE_OPS **ops);

#endif
