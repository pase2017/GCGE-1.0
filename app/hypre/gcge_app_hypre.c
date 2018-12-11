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

void GCGE_HYPRE_VecCreateByVec(void *s_vec, void **d_vec)
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
void GCGE_HYPRE_VecCreateByMat(void *mat, void **vec)
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
void GCGE_HYPRE_VecDestroy(void **vec)
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

void GCGE_HYPRE_LocalVecInnerProd(void *x, void *y, GCGE_DOUBLE *xTy)
{
   hypre_Vector *x_local = hypre_ParVectorLocalVector((HYPRE_ParVector)x);
   hypre_Vector *y_local = hypre_ParVectorLocalVector((HYPRE_ParVector)y);
   *xTy = hypre_SeqVectorInnerProd(x_local, y_local);
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
    ops->VecCreateByVec     = GCGE_HYPRE_VecCreateByVec;
    ops->VecCreateByMat     = GCGE_HYPRE_VecCreateByMat;
    ops->VecDestroy           = GCGE_HYPRE_VecDestroy;

    ops->VecSetRandomValue = GCGE_HYPRE_VecSetRandomValue;
    ops->MatDotVec         = GCGE_HYPRE_MatDotVec;
    ops->VecAxpby          = GCGE_HYPRE_VecAxpby;
    ops->VecInnerProd      = GCGE_HYPRE_VecInnerProd;
    ops->VecLocalInnerProd = GCGE_HYPRE_LocalVecInnerProd;
//    ops->VecLocalInnerProd = GCGE_HYPRE_VecInnerProd;
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
    GCGE_SOLVER_Create(&hypre_solver);
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
        GCGE_HYPRE_VecCreateByMat((void *)A, (void **)evec+i);
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
//    printf("Set up finish!\n");
    return hypre_solver;
}


/**
 * @brief 设置hypre的线性解法器，这里由于hypre没能对解法器进行一个封装，所以必须一个个写，
 *        需要在ops中加入一个成员函数，或者在GCGE_HYPRE_Solver_Init_Linear SolverGivenByUser中新增一个形式参数，
 *        表征用了那个线性解法器，比如，id == 0 表示pcg，id == 1 表示amg, id == 2 表示gemrs
 *        而且对于precond也的setup也会有问题
 *        需要讨论
 *        我认为最佳方案是在GCGE_HYPRE_Solver_Init_Linear SolverGivenByUser中加入新的参数
 *
 * @param Matrix
 * @param b
 * @param x
 * @param ops
 */
void GCGE_HYPRE_GMRESSolve(void *Matrix, void *b, void *x, struct GCGE_OPS_ *ops)
{
   HYPRE_ParCSRFlexGMRESSolve((HYPRE_Solver)(ops->linear_solver_workspace), (HYPRE_ParCSRMatrix)(Matrix), (HYPRE_ParVector)b, (HYPRE_ParVector)x);
}
void GCGE_HYPRE_PCGSolve(void *Matrix, void *b, void *x, struct GCGE_OPS_ *ops)
{
   HYPRE_ParCSRPCGSolve((HYPRE_Solver)(ops->linear_solver_workspace), (HYPRE_ParCSRMatrix)(Matrix), (HYPRE_ParVector)b, (HYPRE_ParVector)x);
}
void GCGE_HYPRE_AMGSolve(void *Matrix, void *b, void *x, struct GCGE_OPS_ *ops)
{
   HYPRE_BoomerAMGSolve((HYPRE_Solver)(ops->linear_solver_workspace), (HYPRE_ParCSRMatrix)(Matrix), (HYPRE_ParVector)b, (HYPRE_ParVector)x);
}

/**
 * @brief 
 *
 * @param A
 * @param B
 * @param linear_solver
 * @param ls_id               表示所用的线性解法器编号，0 AMG, 1 PCG, 2 GMRES
 * @param num_eigenvalues
 * @param argc
 * @param argv[]
 *
 * @return 
 */
GCGE_SOLVER* GCGE_HYPRE_Solver_Init_LinearSolverGivenByUser(HYPRE_ParCSRMatrix A, HYPRE_ParCSRMatrix B, 
      HYPRE_Solver linear_solver, int ls_id, int num_eigenvalues, int argc, char* argv[])
{
    //首先用矩阵A,B创建hypre_solver
    GCGE_SOLVER *hypre_solver = GCGE_HYPRE_Solver_Init(A, B, num_eigenvalues, argc, argv);
    //用户已创建好ksp,直接赋给hypre_solver
    /* 这里也可以直接访问成员函数 */
    GCGE_SOLVER_SetOpsLinearSolverWorkspace(hypre_solver, (void*)linear_solver);
    //将线性解法器设为KSPSolve
    if (ls_id == 0)
    {
//       printf ( "linear solver is AMG\n" );
       hypre_solver->ops->LinearSolver = GCGE_HYPRE_AMGSolve;
    }
    else if (ls_id == 1)
    {
//       printf ( "linear solver is PCG\n" );
       hypre_solver->ops->LinearSolver = GCGE_HYPRE_PCGSolve;
    }
    else if (ls_id == 2)
    {
//       printf ( "linear solver is GMRES\n" );
       hypre_solver->ops->LinearSolver = GCGE_HYPRE_GMRESSolve;
    }
    else {
//       printf ( "Use default CG for linear solver\n" );
    }
    return hypre_solver;
}
