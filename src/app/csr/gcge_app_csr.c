/*
 * =====================================================================================
 *
 *       Filename:  gcge_app_csr.c
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

/* external head file */
#include "gcge_app_csr.h"

void GCGE_CSR_BuildVecByVec(void *s_vec, void **d_vec)
{
   CSR_BuildVecByVec((CSR_VEC*)s_vec, (CSR_VEC**)d_vec);
}
void GCGE_CSR_BuildVecByMat(void *mat, void **vec)
{
   CSR_BuildVecByMat((CSR_MAT*)mat, (CSR_VEC**)vec);
}
void GCGE_CSR_VecFree(void **vec)
{
   CSR_VecFree((CSR_VEC**)vec);
}

void GCGE_CSR_VecSetRandomValue(void *vec)
{
   CSR_VecSetRandomValue((CSR_VEC*)vec);
}
void GCGE_CSR_MatDotVec(void *mat, void *x, void *r)
{
    CSR_MatDotVec((CSR_MAT*)mat, (CSR_VEC*)x, (CSR_VEC*)r);
}
void GCGE_CSR_VecAxpby(GCGE_DOUBLE a, void *x, GCGE_DOUBLE b, void *y)
{
   CSR_VecAxpby(a, (CSR_VEC*)x, b, (CSR_VEC*)y);
}
void GCGE_CSR_VecInnerProd(void *x, void *y, GCGE_DOUBLE *xTy)
{
   CSR_VecInnerProd((CSR_VEC*)x, (CSR_VEC*)y, xTy);
}
void GCGE_CSR_PrintMultiVec(void **x, GCGE_INT n)
{
   GCGE_INT i = 0;
   for(i=0; i<n; i++)
   {
    CSR_PrintVec((CSR_VEC*)x[i]); 
   }
}

void GCGE_CSR_SetOps(GCGE_OPS *ops)
{
    /* either-or */
    ops->BuildVecByVec     = GCGE_CSR_BuildVecByVec;
    ops->BuildVecByMat     = GCGE_CSR_BuildVecByMat;
    ops->FreeVec           = GCGE_CSR_VecFree;

    ops->VecSetRandomValue = GCGE_CSR_VecSetRandomValue;
    ops->MatDotVec         = GCGE_CSR_MatDotVec;
    ops->VecAxpby          = GCGE_CSR_VecAxpby;
    ops->VecInnerProd      = GCGE_CSR_VecInnerProd;
    ops->VecLocalInnerProd = GCGE_CSR_VecInnerProd;
    ops->PrintMultiVec     = GCGE_CSR_PrintMultiVec;
}

void GCGE_SOLVER_SetCSROps(GCGE_SOLVER *solver)
{
    GCGE_CSR_SetOps(solver->ops);
}
