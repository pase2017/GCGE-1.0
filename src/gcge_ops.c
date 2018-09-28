/*
 * =====================================================================================
 *
 *       Filename:  gcge_ops.c
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

#include "gcge_ops.h"


//MultiVec默认值的设定是要依赖Get与Restore这两个函数的
void GCGE_Default_GetVecFromMultiVec(void **V, GCGE_INT j, void **x)
{
    *x = V[j];
}

void GCGE_Default_RestoreVecForMultiVec(void **V, GCGE_INT j, void **x)
{
    V[j] = *x;
}

void GCGE_Default_BuildMultiVecByVec(void *vec, void ***multi_vec, GCGE_INT n_vec, GCGE_OPS *ops)
{
    GCGE_INT i = 0;
    *multi_vec = (void**)malloc(n_vec*sizeof(void*));
    for(i=0; i<n_vec; i++)
    {
        ops->BuildVecByVec(vec, (*multi_vec)+i);
    }
}

void GCGE_Default_BuildMultiVecByMat(void *mat, void ***multi_vec, GCGE_INT n_vec, GCGE_OPS *ops)
{
    GCGE_INT i = 0;
    *multi_vec = (void**)malloc(n_vec*sizeof(void*));
    for(i=0; i<n_vec; i++)
    {
        ops->BuildVecByMat(mat, (*multi_vec)+i);
    }
}

void GCGE_Default_BuildMultiVecByMultiVec(void **init_vec, void ***multi_vec, GCGE_INT n_vec, GCGE_OPS *ops)
{
    GCGE_INT i = 0;
    *multi_vec = (void**)malloc(n_vec*sizeof(void*));
    for(i=0; i<n_vec; i++)
    {
        ops->BuildVecByVec(init_vec[0], (*multi_vec)+i);
    }
}

void GCGE_Default_FreeMultiVec(void ***MultiVec, GCGE_INT n_vec, struct GCGE_OPS_ *ops)
{
    GCGE_INT i = 0;
    for(i=0; i<n_vec; i++)
    {
        ops->FreeVec(&(*MultiVec)[i]);  (*MultiVec)[i] = NULL;
    }
    free(*MultiVec); *MultiVec = NULL;
}

void GCGE_Default_MultiVecSetRandomValue(void **multi_vec, GCGE_INT n_vec, struct GCGE_OPS_ *ops)
{
    GCGE_INT i = 0;
    void     *x;
    for(i=0; i<n_vec; i++)
    {
        ops->GetVecFromMultiVec(multi_vec, i, &x);
        ops->VecSetRandomValue(x);
        ops->RestoreVecForMultiVec(multi_vec, i, &x);
    }
}

void GCGE_Default_MatDotMultiVec(void *mat, void **x, void **y, 
        GCGE_INT *start, GCGE_INT *end, struct GCGE_OPS_ *ops)
{
    GCGE_INT i = 0;
    GCGE_INT n_vec = end[0]-start[0];
    void     *xs, *ys;
    for(i=0; i<n_vec; i++)
    {
        ops->GetVecFromMultiVec(x, start[0]+i, &xs);
        ops->GetVecFromMultiVec(y, start[1]+i, &ys);
        ops->MatDotVec(mat, x[start[0]+i], y[start[1]+i]);
        ops->RestoreVecForMultiVec(x, start[0]+i, &xs);
        ops->RestoreVecForMultiVec(y, start[1]+i, &ys);
    }
}

void GCGE_Default_MultiVecAxpby(GCGE_DOUBLE a, void **x, GCGE_DOUBLE b, void **y, 
        GCGE_INT *start, GCGE_INT *end, struct GCGE_OPS_ *ops)
{
    GCGE_INT i = 0;
    GCGE_INT n_vec = end[0]-start[0];
    void     *xs, *ys;
    for(i=0; i<n_vec; i++)
    {
        ops->GetVecFromMultiVec(x, start[0]+i, &xs);
        ops->GetVecFromMultiVec(y, start[1]+i, &ys);
        ops->VecAxpby(a, x[start[0]+i], b, y[start[1]+i]);
        ops->RestoreVecForMultiVec(x, start[0]+i, &xs);
        ops->RestoreVecForMultiVec(y, start[1]+i, &ys);
    }
}

void GCGE_Default_MultiVecAxpbyColumn(GCGE_DOUBLE a, void **x, GCGE_INT col_x, 
        GCGE_DOUBLE b, void **y, GCGE_INT col_y, struct GCGE_OPS_ *ops)
{
    void     *xs, *ys;
    ops->GetVecFromMultiVec(x, col_x, &xs);
    ops->GetVecFromMultiVec(y, col_y, &ys);
    ops->VecAxpby(a, x[col_x], b, y[col_y]);
    ops->RestoreVecForMultiVec(x, col_x, &xs);
    ops->RestoreVecForMultiVec(y, col_y, &ys);
}

/* vec_y[j] = \sum_{i=sx}^{ex} vec_x[i] a[i][j] */
void GCGE_Default_MultiVecLinearComb(void **x, void **y, GCGE_INT *start, GCGE_INT *end,
        GCGE_DOUBLE *a, GCGE_INT lda, void *dmat, GCGE_INT lddmat, struct GCGE_OPS_ *ops)
{
    GCGE_INT idx_x = 0;
    GCGE_INT idx_y = 0;
    GCGE_INT nx = end[0]-start[0];
    GCGE_INT ny = end[1]-start[1];
    void     *xs, *ys;
    for(idx_y=start[1]; idx_y<end[1]; idx_y++)
    {
        ops->GetVecFromMultiVec(x, start[0], &xs);
        ops->GetVecFromMultiVec(y, idx_y,    &ys);

        ops->VecAxpby(a[idx_y*lda+start[0]], xs, 0.0, ys);

        ops->RestoreVecForMultiVec(x, start[0], &xs);
        
        for(idx_x = start[0]+1; idx_x<end[0]; idx_x++)
        {
            ops->GetVecFromMultiVec(x, idx_x, &xs);
            ops->VecAxpby(a[idx_y*lda+idx_x], xs, 1.0, ys);
            ops->RestoreVecForMultiVec(x, idx_x, &xs);
        }
        ops->RestoreVecForMultiVec(y, idx_y, &ys);
    }
}

void GCGE_Default_MultiVecInnerProd(void **V, void **W, GCGE_DOUBLE *a, char *is_sym, 
        GCGE_INT *start, GCGE_INT *end, GCGE_INT lda, struct GCGE_OPS_ *ops)
{
    GCGE_INT iv = 0, iw = 0;
    GCGE_DOUBLE *p = a;
    void     *vs, *ws;
    /* TODO some error */
    if(strcmp(is_sym, "sym") == 0)
    {
        GCGE_INT length = end[0] - start[0];
        for(iw=0, p=a; iw<length; iw++)
        {
            ops->GetVecFromMultiVec(W, iw+start[1], &ws);
            for(iv=0; iv<=iw; iv++, p++)
            {
                ops->GetVecFromMultiVec(V, iv+start[0], &vs);
                ops->VecInnerProd(vs, ws, p);
                ops->RestoreVecForMultiVec(V, iv+start[0], &vs);
//		printf ( "iw = %d, iv = %d\n", iw, iv );
            }
            ops->RestoreVecForMultiVec(W, iw+start[1], &ws);
            p += lda-(iw+1);
        }
        //对称化
//	printf ( "-----------------\n" );
        for(iw=0; iw<length; iw++)
        {
            for(iv=iw+1; iv<length; iv++)
            {
//		printf ( "iw = %d, iv = %d\n", iw, iv );
                a[iw*lda+iv] = a[iv*lda+iw];
            }
        }
    }
    else
    {
        GCGE_INT jvw = lda - (end[0] - start[0]);
        for(iw=start[1], p=a; iw<end[1]; iw++)
        {
            ops->GetVecFromMultiVec(W, iw, &ws);
            for(iv=start[0]; iv<end[0]; iv++, p++)
            {
                ops->GetVecFromMultiVec(V, iv, &vs);
                ops->VecInnerProd(V[iv], W[iw], p);
                ops->RestoreVecForMultiVec(V, iv, &vs);
            }
            ops->RestoreVecForMultiVec(W, iw, &ws);
            p += jvw;
        }
    }

}

void GCGE_Default_MultiVecSwap(void **V_1, void **V_2, GCGE_INT *start, GCGE_INT *end, 
        struct GCGE_OPS_ *ops)
{
    GCGE_INT i = 0;
    GCGE_INT size = end[0]-start[0];
    void     *vs, *ws;
    for(i=0; i<size; i++)
    {
        ops->GetVecFromMultiVec(V_1, start[0]+i, &vs);
        ops->GetVecFromMultiVec(V_2, start[1]+i, &ws);
        ops->RestoreVecForMultiVec(V_1, start[0]+i, &ws);
        ops->RestoreVecForMultiVec(V_2, start[1]+i, &vs);
    }
}

void GCGE_Default_DenseMatEigenSolver(char *jobz, char *range, char *uplo, 
        GCGE_INT *nrows, GCGE_DOUBLE *a, GCGE_INT *lda, 
        GCGE_DOUBLE *vl, GCGE_DOUBLE *vu, GCGE_INT *il, GCGE_INT *iu, 
        GCGE_DOUBLE *abstol, GCGE_INT *nev, 
        GCGE_DOUBLE *eval, GCGE_DOUBLE *evec, GCGE_INT *lde, 
        GCGE_INT *isuppz, GCGE_DOUBLE *work, GCGE_INT *lwork, 
        GCGE_INT *iwork, GCGE_INT *liwork, 
        GCGE_INT *ifail, GCGE_INT *info)
{
    dsyevx_(jobz, range, uplo, nrows, a, lda, vl, vu, il, iu, 
            abstol, nev, eval, evec, lde, work, lwork, 
            iwork, ifail, info);
    
    /*
    *lwork = 26*nrows;
    *liwork = 10*nrows;
    dsyevr_(jobz, range, uplo, nrows, a, lda, vl, vu, il, iu, 
            abstol, nev, eval, evec, lde, isuppz, work, lwork, 
            iwork, liwork, info);
            */

    /*
    *lwork = 3*nrows;
    dsyev_(jobz, uplo, nrows, a, lda, eval, work, lwork, info);
    memcpy(evec, a, nrows*nrows*sizeof(GCGE_DOUBLE));
    */

    /*
    *lwork = 1+6*nrows+2*nrows*nrows;
    dsyevd_(jobz, uplo, nrows, a, lda, eval, work, iwork, info);
    memcpy(evec, a, nrows*nrows*sizeof(GCGE_DOUBLE));
    */
}

void GCGE_Default_DenseMatDotDenseMat(char *transa, char *transb, 
        GCGE_INT    *nrows, GCGE_INT    *ncols, GCGE_INT    *mid, 
        GCGE_DOUBLE *alpha, GCGE_DOUBLE *a,     GCGE_INT    *lda, 
        GCGE_DOUBLE *b,     GCGE_INT    *ldb,   GCGE_DOUBLE *beta, 
        GCGE_DOUBLE *c,     GCGE_INT    *ldc)
{
    dgemm_(transa, transb, nrows, ncols, mid, 
            alpha, a, lda, b, ldb, beta, c, ldc);
}

void GCGE_Default_DenseSymMatDotDenseMat(char *side, char *uplo, 
        GCGE_INT    *nrows, GCGE_INT    *ncols,
        GCGE_DOUBLE *alpha, GCGE_DOUBLE *a,   GCGE_INT    *lda,
        GCGE_DOUBLE *b,     GCGE_INT    *ldb, GCGE_DOUBLE *beta,
        GCGE_DOUBLE *c,     GCGE_INT    *ldc)
{
    dsymm_(side, uplo, nrows, ncols, alpha, a, lda,
           b, ldb, beta, c, ldc);
}


void GCGE_OPS_Create(GCGE_OPS **ops)
{
    *ops = (GCGE_OPS*)malloc(sizeof(GCGE_OPS));
    (*ops)->VecSetRandomValue = NULL;
    (*ops)->MatDotVec = NULL;
    (*ops)->VecAxpby = NULL;
    (*ops)->VecInnerProd = NULL;
    (*ops)->VecLocalInnerProd = NULL;
    (*ops)->BuildVecByVec = NULL;
    (*ops)->FreeVec = NULL;
    (*ops)->BuildMultiVecByVec = NULL;
    (*ops)->BuildMultiVecByMat = NULL;
    (*ops)->BuildMultiVecByMultiVec = NULL;
    (*ops)->FreeMultiVec = NULL;
    (*ops)->MultiVecSetRandomValue = NULL;
    (*ops)->MatDotMultiVec = NULL;
    (*ops)->MultiVecAxpby = NULL;
    (*ops)->MultiVecAxpbyColumn = NULL;
    (*ops)->MultiVecLinearComb = NULL;
    (*ops)->MultiVecInnerProd = NULL;
    (*ops)->MultiVecSwap = NULL;
    (*ops)->GetVecFromMultiVec = NULL;
    (*ops)->RestoreVecForMultiVec = NULL;
    (*ops)->DenseMatEigenSolver = NULL;
    (*ops)->DenseMatDotDenseMat = NULL;
    (*ops)->DenseSymMatDotDenseMat = NULL;
    (*ops)->Orthogonal = NULL;
    (*ops)->DenseMatCreate = NULL;
    (*ops)->DenseMatDestroy = NULL;
    (*ops)->LinearSolver = NULL;
}

GCGE_INT GCGE_OPS_Setup(GCGE_OPS *ops)
{
    GCGE_INT error;
    if(ops->MatDotVec == NULL)
    {
        error = 1;
        printf("ERROR: Please provide the MatDotVec operation!\n");
        return error;
    }
    if(ops->VecAxpby == NULL)
    {
        error = 1;
        printf("ERROR: Please provide the VecAxpby operation!\n");
        return error;
    }
    if(ops->VecInnerProd == NULL)
    {
        error = 1;
        printf("ERROR: Please provide the VecInnerProd operation!\n");
        return error;
    }
    if(ops->VecLocalInnerProd == NULL)
    {
        error = 1;
        printf("ERROR: Please provide the VecLocalInnerProd operation!\n");
        return error;
    }
    if(ops->BuildVecByVec == NULL)
    {
        error = 1;
        printf("ERROR: Please provide the BuildVecByVec operation!\n");
        return error;
    }
    if(ops->BuildVecByMat == NULL)
    {
        error = 1;
        printf("ERROR: Please provide the BuildVecByMat operation!\n");
        return error;
    }
    if(ops->FreeVec == NULL)
    {
        error = 1;
        printf("ERROR: Please provide the FreeVec operation!\n");
        return error;
    }
    if(ops->VecSetRandomValue == NULL)
    {
        error = 1;
        printf("ERROR: Please provide the VecSetRandomValue operation!\n");
        return error;
    }
    if(ops->BuildMultiVecByVec == NULL)
    {
        ops->BuildMultiVecByVec = GCGE_Default_BuildMultiVecByVec;
    }
    if(ops->BuildMultiVecByMat == NULL)
    {
        ops->BuildMultiVecByMat = GCGE_Default_BuildMultiVecByMat;
    }
    if(ops->BuildMultiVecByMultiVec == NULL)
    {
        ops->BuildMultiVecByMultiVec = GCGE_Default_BuildMultiVecByMultiVec;
    }
    if(ops->FreeMultiVec == NULL)
    {
        ops->FreeMultiVec = GCGE_Default_FreeMultiVec; 
    }
    if(ops->MultiVecSetRandomValue == NULL)
    {
        ops->MultiVecSetRandomValue = GCGE_Default_MultiVecSetRandomValue;
    }
    if(ops->MatDotMultiVec == NULL)
    {
        ops->MatDotMultiVec = GCGE_Default_MatDotMultiVec;
    }
    if(ops->MultiVecAxpby == NULL)
    {
        ops->MultiVecAxpby = GCGE_Default_MultiVecAxpby;
    }
    if(ops->MultiVecAxpbyColumn == NULL)
    {
        ops->MultiVecAxpbyColumn = GCGE_Default_MultiVecAxpbyColumn;
    }
    if(ops->MultiVecLinearComb == NULL)
    {
        ops->MultiVecLinearComb = GCGE_Default_MultiVecLinearComb;
    }
    if(ops->MultiVecInnerProd == NULL)
    {
        ops->MultiVecInnerProd = GCGE_Default_MultiVecInnerProd;
    }
    if(ops->MultiVecSwap == NULL)
    {
        ops->MultiVecSwap = GCGE_Default_MultiVecSwap;
    }
    if(ops->GetVecFromMultiVec == NULL)
    {
        ops->GetVecFromMultiVec = GCGE_Default_GetVecFromMultiVec;
    }
    if(ops->RestoreVecForMultiVec == NULL)
    {
        ops->RestoreVecForMultiVec = GCGE_Default_RestoreVecForMultiVec;
    }
    if(ops->DenseMatEigenSolver == NULL)
    {
        ops->DenseMatEigenSolver = GCGE_Default_DenseMatEigenSolver;
    }
    if(ops->DenseMatDotDenseMat == NULL)
    {
        ops->DenseMatDotDenseMat = GCGE_Default_DenseMatDotDenseMat;
    }
    if(ops->DenseSymMatDotDenseMat == NULL)
    {
        ops->DenseSymMatDotDenseMat = GCGE_Default_DenseSymMatDotDenseMat;
    }
}                                           
                                            
void GCGE_OPS_Free(GCGE_OPS **ops)
{
    free(*ops); *ops = NULL;
}


