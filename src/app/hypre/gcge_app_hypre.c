/*
 * =====================================================================================
 *
 *       Filename:  gcge_app_hypre.c
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
#include "gcge_app_hypre.h"

void GCGE_HYPRE_BuildVecByVec(void *s_vec, void **d_vec)
{
  HYPRE_ParVector  x_hypre      = (HYPRE_ParVector)s_vec;
  MPI_Comm         comm         = hypre_ParVectorComm(x_hypre);
  GCGE_INT         global_size  = hypre_ParVectorGlobalSize(x_hypre);
  GCGE_INT        *partitioning = hypre_ParVectorPartitioning(x_hypre);
//  printf ( "vec by vec partitioning: %d, %d, global_size = %d\n", partitioning[0], partitioning[1], global_size );
  HYPRE_ParVector  y_hypre      = hypre_ParVectorCreate(comm, global_size, partitioning);
  hypre_ParVectorInitialize(y_hypre);
  hypre_ParVectorSetPartitioningOwner(y_hypre, 0);
  *d_vec = (void *)y_hypre;
}
void GCGE_HYPRE_BuildVecByMat(void *mat, void **vec)
{
  HYPRE_ParCSRMatrix A_hypre      = (HYPRE_ParCSRMatrix)mat;
  MPI_Comm           comm         = hypre_ParCSRMatrixComm(A_hypre);
  GCGE_INT           global_size  = hypre_ParCSRMatrixGlobalNumRows(A_hypre);
  GCGE_INT          *partitioning = NULL;
#ifdef HYPRE_NO_GLOBAL_PARTITION
  partitioning = hypre_CTAlloc(HYPRE_Int,  2);
  hypre_ParCSRMatrixGetLocalRange(A_hypre, partitioning, partitioning+1, partitioning, partitioning+1);
  partitioning[1] += 1;
#else
  HYPRE_ParCSRMatrixGetRowPartitioning(A_hypre, &partitioning);
#endif
//  printf ( "vec by mat partitioning: %d, %d, global_size = %d\n", partitioning[0], partitioning[1], global_size );
  HYPRE_ParVector    y_hypre      = hypre_ParVectorCreate(comm, global_size, partitioning);
  hypre_ParVectorInitialize(y_hypre);
  hypre_ParVectorSetPartitioningOwner(y_hypre, 1);
  *vec = (void *)y_hypre;
}
void GCGE_HYPRE_VecFree(void **vec)
{
   hypre_ParVectorDestroy((HYPRE_ParVector)(*vec));
   *vec = NULL;
}

void GCGE_HYPRE_VecSetRandomValue(void *vec)
{
   hypre_ParVectorSetRandomValues((HYPRE_ParVector)vec, rand());
}
void GCGE_HYPRE_MatDotVec(void *mat, void *x, void *r)
{
   hypre_ParCSRMatrixMatvec(1.0,  (HYPRE_ParCSRMatrix)mat,  (HYPRE_ParVector)x,  0.0,  (HYPRE_ParVector)r);
}
void GCGE_HYPRE_VecAxpby(GCGE_DOUBLE a, void *x, GCGE_DOUBLE b, void *y)
{
   hypre_ParVectorScale(b, (HYPRE_ParVector)y);
   hypre_ParVectorAxpy(a, (HYPRE_ParVector)x, (HYPRE_ParVector)y);
}
void GCGE_HYPRE_VecInnerProd(void *x, void *y, GCGE_DOUBLE *xTy)
{
   HYPRE_ParVectorInnerProd((HYPRE_ParVector)x,  (HYPRE_ParVector)y,  xTy);
}
void GCGE_HYPRE_PrintMultiVec(void **x, GCGE_INT n)
{
   char file_name[16];
   GCGE_INT i = 0;
   for(i=0; i<n; i++)
   {
      hypre_sprintf(file_name, "%s.%d", "multivec", i);
      hypre_ParVectorPrint((HYPRE_ParVector)x[i],  file_name);
   }
}

void GCGE_HYPRE_SetOps(GCGE_OPS *ops)
{
    /* either-or */
    ops->BuildVecByVec     = GCGE_HYPRE_BuildVecByVec;
    ops->BuildVecByMat     = GCGE_HYPRE_BuildVecByMat;
    ops->FreeVec           = GCGE_HYPRE_VecFree;

    ops->VecSetRandomValue = GCGE_HYPRE_VecSetRandomValue;
    ops->MatDotVec         = GCGE_HYPRE_MatDotVec;
    ops->VecAxpby          = GCGE_HYPRE_VecAxpby;
    ops->VecInnerProd      = GCGE_HYPRE_VecInnerProd;
    ops->VecLocalInnerProd = GCGE_HYPRE_VecInnerProd;
    ops->PrintMultiVec     = GCGE_HYPRE_PrintMultiVec;
}

void GCGE_SOLVER_SetHYPREOps(GCGE_SOLVER *solver)
{
    GCGE_HYPRE_SetOps(solver->ops);
}

//下面是一个对HYPRE 的GCG_Solver的初始化
GCGE_SOLVER* GCGE_HYPRE_Solver_Init(HYPRE_ParCSRMatrix A, HYPRE_ParCSRMatrix B, int num_eigenvalues, int argc, char* argv[])
{
    //第一步: 定义相应的ops,para以及workspace

    //创建一个solver变量
    GCGE_SOLVER *hypre_solver;
    GCGE_INT help = GCGE_SOLVER_Create(&hypre_solver);
    if(help)
    {
        GCGE_SOLVER_Free(&hypre_solver);
        exit(0);
    }
    if(num_eigenvalues != -1)
        hypre_solver->para->nev = num_eigenvalues;
    //设置初始值
    GCGE_INT error = GCGE_PARA_SetFromCommandLine(hypre_solver->para, argc, argv);
    int nev = hypre_solver->para->nev;
    double *eval = (double *)calloc(nev, sizeof(double)); 
    hypre_solver->eval = eval;

    HYPRE_ParVector *evec = hypre_CTAlloc(HYPRE_ParVector,  nev);
    int i = 0; 
    //这里为什么不一次性生成 CSR_BuildMultiVecByMat?
    for(i=0; i<nev; i++)
    {
        GCGE_HYPRE_BuildVecByMat((void *)A, (void **)evec+i);
    }
    hypre_solver->evec = (void**)evec;

    GCGE_SOLVER_SetMatA(hypre_solver, A);
    if(B != NULL)
        GCGE_SOLVER_SetMatB(hypre_solver, B);
    GCGE_SOLVER_SetHYPREOps(hypre_solver);
    GCGE_SOLVER_SetEigenvalues(hypre_solver, eval);
    GCGE_SOLVER_SetEigenvectors(hypre_solver, (void**)evec);
    //setup and solve
    GCGE_SOLVER_Setup(hypre_solver);
    printf("Set up finish!\n");
    return hypre_solver;
}
