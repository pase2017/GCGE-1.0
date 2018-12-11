/*
 * =====================================================================================
 *
 *       Filename:  gcge_app_phg.c
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
#include "gcge_app_phg.h"

void GCGE_PHG_VecCreateByVec(void *s_vec, void **d_vec)
{
   VEC* x_phg = (VEC*)s_vec;
   VEC* y_phg = phgMapCreateVec(x_phg->map, 1);
   *d_vec = (void *)y_phg;
}
void GCGE_PHG_VecCreateByMat(void *mat, void **vec)
{
  MAT* A_phg = (MAT*)mat;
  VEC* y_phg = phgMapCreateVec(A_phg->cmap, 1);
  *vec = (void *)y_phg;
}
void GCGE_PHG_VecDestroy(void **vec)
{
   phgVecDestroy((VEC**)(vec));
}

void GCGE_PHG_VecSetRandomValue(void *vec)
{
   phgVecRandomize((VEC *)vec, rand());
}
void GCGE_PHG_MatDotVec(void *mat, void *x, void *r)
{
   phgMatVec(0, 1.0, (MAT *)mat, (VEC *)x, 0.0, (VEC **)&r);
}
void GCGE_PHG_VecAxpby(GCGE_DOUBLE a, void *x, GCGE_DOUBLE b, void *y)
{
   phgMatVec(0, a, NULL, (VEC *)x, b, (VEC **)&y);
}
void GCGE_PHG_VecInnerProd(void *x, void *y, GCGE_DOUBLE *xTy)
{
   phgVecDot((VEC *)x, 0, (VEC *)y, 0,  xTy);
}

void GCGE_PHG_LocalVecInnerProd(void *x, void *y, GCGE_DOUBLE *xTy)
{
   int i;
   *xTy = 0;
   for (i = 0; i < ((VEC *)x)->map->nlocal; i++)
      *xTy += ((VEC *)x)->data[i] * ((VEC *)y)->data[i];
}
void GCGE_PHG_MultiVecPrint(void **x, GCGE_INT n)
{
   char file_name[16];
   GCGE_INT i = 0;
   for(i=0; i<n; i++)
   {
      phgVecDump((VEC *)x[i],  "multivec");
   }
}
void GCGE_PHG_SetOps(GCGE_OPS *ops)
{
    /* either-or */
    ops->VecCreateByVec     = GCGE_PHG_VecCreateByVec;
    ops->VecCreateByMat     = GCGE_PHG_VecCreateByMat;
    ops->VecDestroy           = GCGE_PHG_VecDestroy;

    ops->VecSetRandomValue = GCGE_PHG_VecSetRandomValue;
    ops->MatDotVec         = GCGE_PHG_MatDotVec;
    ops->VecAxpby          = GCGE_PHG_VecAxpby;
    ops->VecInnerProd      = GCGE_PHG_VecInnerProd;
    ops->VecLocalInnerProd = GCGE_PHG_LocalVecInnerProd;
//    ops->VecLocalInnerProd = GCGE_PHG_VecInnerProd;
    ops->MultiVecPrint     = GCGE_PHG_MultiVecPrint;
}

void GCGE_SOLVER_SetPHGOps(GCGE_SOLVER *solver)
{
    GCGE_PHG_SetOps(solver->ops);
}

//下面是一个对PHG 的GCG_Solver的初始化
GCGE_SOLVER* GCGE_PHG_Solver_Init(MAT *A, MAT *B, int num_eigenvalues, int argc, char* argv[])
{
    //第一步: 定义相应的ops,para以及workspace

    //创建一个solver变量
    GCGE_SOLVER *phg_solver;
    GCGE_SOLVER_Create(&phg_solver);
    if(num_eigenvalues != -1)
        phg_solver->para->nev = num_eigenvalues;
    //设置初始值
    GCGE_INT error = GCGE_PARA_SetFromCommandLine(phg_solver->para, argc, argv);
    int nev = phg_solver->para->nev;
    double *eval = (double *)calloc(nev, sizeof(double)); 
    phg_solver->eval = eval;

    VEC** evec = phgCalloc(sizeof(VEC*),  nev);
    int i = 0; 
    //这里为什么不一次性生成 CSR_MultiVecCreateByMat?
    for(i=0; i<nev; i++)
    {
        GCGE_PHG_VecCreateByMat((void *)A, (void **)evec+i);
    }
    phg_solver->evec = (void**)evec;

    GCGE_SOLVER_SetMatA(phg_solver, A);
    if(B != NULL)
        GCGE_SOLVER_SetMatB(phg_solver, B);
    GCGE_SOLVER_SetPHGOps(phg_solver);
    GCGE_SOLVER_SetEigenvalues(phg_solver, eval);
    GCGE_SOLVER_SetEigenvectors(phg_solver, (void**)evec);
    //setup and solve
    GCGE_SOLVER_Setup(phg_solver);
//    printf("Set up finish!\n");
    return phg_solver;
}


/**
 * @brief 设置phg的线性解法器，
 *
 * @param Matrix
 * @param b
 * @param x
 * @param ops
 */
void GCGE_PHG_Solve(void *Matrix, void *b, void *x, struct GCGE_OPS_ *ops)
{
   SOLVER* solver = ops->linear_solver_workspace;

   solver->rhs->data = ((VEC *)b)->data;
   solver->rhs->assembled = TRUE;
   phgSolverVecSolve(solver,  FALSE, (VEC *)x);
   solver->rhs->data = NULL;
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
GCGE_SOLVER* GCGE_PHG_Solver_Init_LinearSolverGivenByUser(MAT *A, MAT *B, 
      SOLVER* linear_solver, int num_eigenvalues, int argc, char* argv[])
{
    //首先用矩阵A,B创建phg_solver
    GCGE_SOLVER *phg_solver = GCGE_PHG_Solver_Init(A, B, num_eigenvalues, argc, argv);
    //用户已创建好ksp,直接赋给phg_solver
    /* 这里也可以直接访问成员函数 */
    GCGE_SOLVER_SetOpsLinearSolverWorkspace(phg_solver, (void*)linear_solver);
    phg_solver->ops->LinearSolver = GCGE_PHG_Solve;

    return phg_solver;
}
