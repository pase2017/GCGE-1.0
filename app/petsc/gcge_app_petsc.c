/*
 * =====================================================================================
 *
 *       Filename:  gcge_app_petsc.c
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
#include <time.h>

#include "gcge_app_petsc.h"

void PETSC_ReadMatrixBinary(Mat *A, const char *filename)
{
    PetscErrorCode ierr;
    PetscViewer    viewer;
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Getting matrix...\n"); 
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer); 
    ierr = MatCreate(PETSC_COMM_WORLD, A); 
    ierr = MatSetFromOptions(*A); 
    ierr = MatLoad(*A, viewer); 
    ierr = PetscViewerDestroy(&viewer);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Getting matrix... Done\n");
}

void PETSC_LinearSolverCreate(KSP *ksp, Mat A, Mat T)
{
    PetscErrorCode ierr;
    ierr = KSPCreate(PETSC_COMM_WORLD,ksp);
    ierr = KSPSetOperators(*ksp,A,T);
    ierr = KSPSetInitialGuessNonzero(*ksp, 1);

    ierr = KSPSetType(*ksp, KSPCG);
    //这里的rtol应取作<=ev_tol
    //PetscErrorCode  KSPSetTolerances(KSP ksp,PetscReal rtol,PetscReal abstol,PetscReal dtol,PetscInt maxits)
    ierr = KSPSetTolerances(*ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, 1000);
    PC pc;
    ierr = KSPGetPC(*ksp, &pc);
    ierr = PCSetType(pc, PCHYPRE);
    ierr = PCHYPRESetType(pc, "boomeramg");
    //最后从命令行设置参数
    ierr = KSPSetFromOptions(*ksp);
}

void PETSC_VecLocalInnerProd(Vec x, Vec y, double *value)
{
    PetscErrorCode     ierr;
    const PetscScalar *local_x;
    const PetscScalar *local_y;
    PetscInt           low, high, length, i = 0;
    ierr = VecGetOwnershipRange(x, &low, &high);
    length = high-low;
    ierr = VecGetArrayRead(x,&local_x);
    ierr = VecGetArrayRead(y,&local_y);
    *value = 0.0;
    for(i=0; i<length; i++)
    {
        *value += local_x[i]*local_y[i];
    }
    ierr = VecRestoreArrayRead(x,&local_x);
    ierr = VecRestoreArrayRead(y,&local_y);
}


void GCGE_PETSC_VecCreateByVec(void **d_vec, void *s_vec)
{
	PetscErrorCode ierr;
    ierr = VecDuplicate((Vec)s_vec, (Vec*)d_vec);
}
void GCGE_PETSC_VecCreateByMat(void **vec, void *mat)
{
	PetscErrorCode ierr;
    ierr = MatCreateVecs((Mat)mat, NULL, (Vec*)vec);
}
void GCGE_PETSC_MultiVecCreateByMat(void ***multi_vec, GCGE_INT n_vec, void *mat, GCGE_OPS *ops)
{
	PetscErrorCode ierr;
	Vec vector;
    ierr = MatCreateVecs((Mat)mat, NULL, &vector);
    ierr = VecDuplicateVecs(vector, n_vec, (Vec**)multi_vec);
}
void GCGE_PETSC_VecDestroy(void **vec)
{
	PetscErrorCode ierr;
    ierr = VecDestroy((Vec*)vec);
}
void GCGE_PETSC_MultiVecDestroy(void ***MultiVec, GCGE_INT n_vec, struct GCGE_OPS_ *ops)
{
	PetscErrorCode ierr;
    ierr = VecDestroyVecs(n_vec, (Vec**)MultiVec);
}

void GCGE_PETSC_VecSetRandomValue(void *vec)
{
    PetscRandom    rctx;
    PetscErrorCode ierr;
    ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
    ierr = PetscRandomSetFromOptions(rctx);
    ierr = VecSetRandom((Vec)vec, rctx);
    ierr = PetscRandomDestroy(&rctx);
}
void GCGE_PETSC_MultiVecSetRandomValue(void **multi_vec, GCGE_INT start, GCGE_INT n_vec, struct GCGE_OPS_ *ops)
{
    PetscRandom    rctx;
    PetscErrorCode ierr;
    PetscInt       i;
    ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
    ierr = PetscRandomSetFromOptions(rctx);
    for(i=0; i<n_vec; i++)
    {
        ierr = VecSetRandom(((Vec*)multi_vec)[i], rctx);
    }
    ierr = PetscRandomDestroy(&rctx);
}

void GCGE_PETSC_MatDotVec(void *mat, void *x, void *r)
{
    PetscErrorCode ierr = MatMult((Mat)mat, (Vec)x, (Vec)r);
}
void GCGE_PETSC_VecAxpby(GCGE_DOUBLE a, void *x, GCGE_DOUBLE b, void *y)
{
    PetscErrorCode ierr;
	ierr = VecScale((Vec)y, b);
    if(x != y)
    {
	    ierr = VecAXPY((Vec)y, a, (Vec)x);
    }
}
void GCGE_PETSC_VecInnerProd(void *x, void *y, GCGE_DOUBLE *value_ip)
{
	PetscErrorCode ierr = VecDot((Vec)x, (Vec)y, value_ip);
}

void GCGE_PETSC_VecLocalInnerProd(void *x, void *y, GCGE_DOUBLE *value_ip)
{
    PETSC_VecLocalInnerProd((Vec)x, (Vec)y, value_ip);
}

void GCGE_PETSC_LinearSolver(void *Matrix, void *b, void *x, struct GCGE_OPS_ *ops)
{
    PetscErrorCode ierr;
    ierr = KSPSolve((KSP)(ops->linear_solver_workspace), (Vec)b, (Vec)x);
}

void GCGE_PETSC_SetOps(GCGE_OPS *ops)
{
    /* either-or */
    ops->VecCreateByVec     = GCGE_PETSC_VecCreateByVec;
    ops->VecCreateByMat     = GCGE_PETSC_VecCreateByMat;
    ops->VecDestroy           = GCGE_PETSC_VecDestroy;

    ops->VecSetRandomValue = GCGE_PETSC_VecSetRandomValue;
    ops->MatDotVec         = GCGE_PETSC_MatDotVec;
    ops->VecAxpby          = GCGE_PETSC_VecAxpby;
    ops->VecInnerProd      = GCGE_PETSC_VecInnerProd;
    ops->VecLocalInnerProd = GCGE_PETSC_VecLocalInnerProd;

    ops->MultiVecDestroy           = GCGE_PETSC_MultiVecDestroy;
    ops->MultiVecCreateByMat     = GCGE_PETSC_MultiVecCreateByMat;
    ops->MultiVecSetRandomValue = GCGE_PETSC_MultiVecSetRandomValue;
}

void GCGE_SOLVER_SetPETSCOps(GCGE_SOLVER *solver)
{
    GCGE_PETSC_SetOps(solver->ops);
}

//下面是一个对PETSC 的GCG_Solver的初始化
//用户只提供要求解特征值问题的矩阵A,B,线性解法器等都用默认值
//如果用户想在命令行提供特征值个数，需要将num_eigenvalues赋值为-1
GCGE_SOLVER* GCGE_PETSC_Solver_Init(Mat A, Mat B, int num_eigenvalues, int argc, char* argv[])
{
    //第一步: 定义相应的ops,para以及workspace

    //创建一个solver变量
    GCGE_SOLVER *petsc_solver;
    GCGE_SOLVER_Create(&petsc_solver);
    if(num_eigenvalues != -1)
        petsc_solver->para->nev = num_eigenvalues;
    GCGE_INT error = GCGE_PARA_SetFromCommandLine(petsc_solver->para, argc, argv);
    //设置初始值
    int nev = petsc_solver->para->nev;
    double *eval = (double *)calloc(nev, sizeof(double)); 
    petsc_solver->eval = eval;

	PetscErrorCode ierr;
    Vec *evec;
	Vec vector;
    ierr = MatCreateVecs(A, NULL, &vector);
    ierr = VecDuplicateVecs(vector, nev, &evec);

    GCGE_SOLVER_SetMatA(petsc_solver, A);
    if(B != NULL)
        GCGE_SOLVER_SetMatB(petsc_solver, B);
    GCGE_SOLVER_SetPETSCOps(petsc_solver);
    GCGE_SOLVER_SetEigenvalues(petsc_solver, eval);
    GCGE_SOLVER_SetEigenvectors(petsc_solver, (void**)evec);
    //setup and solve
    GCGE_SOLVER_Setup(petsc_solver);
    GCGE_PrintParaInfo(petsc_solver->para);
    GCGE_Printf("Set up finish!\n");
    return petsc_solver;
}

//下面是一个对PETSC 的GCG_Solver的初始化
//用户只提供要求解特征值问题的矩阵A,B,线性解法器使用KSPSolve,用户不接触ksp,设置参数均从命令行设置
GCGE_SOLVER* GCGE_PETSC_Solver_Init_KSPDefault(Mat A, Mat B, Mat P, int num_eigenvalues, int argc, char* argv[])
{
    //首先用矩阵A,B创建petsc_solver
    GCGE_SOLVER *petsc_solver = GCGE_PETSC_Solver_Init(A, B, num_eigenvalues, argc, argv);
    KSP ksp;
    //使用矩阵A和预条件矩阵P创建ksp,ksp的参数均为默认参数,要修改ksp参数只能用命令行参数进行设置
    PETSC_LinearSolverCreate(&ksp, (Mat)A, (Mat)P);
    //把ksp作为ops->linear_solver_workspace
    GCGE_SOLVER_SetOpsLinearSolverWorkspace(petsc_solver, (void*)ksp);
    //将线性解法器设为KSPSolve
    petsc_solver->ops->LinearSolver = GCGE_PETSC_LinearSolver;
    return petsc_solver;
}

//下面是一个对PETSC 的GCG_Solver的初始化
//用户提供要求解特征值问题的矩阵A,B,以及线性解法器结构ksp,使用KSPSolve作为线性解法器
GCGE_SOLVER* GCGE_PETSC_Solver_Init_KSPGivenByUser(Mat A, Mat B, KSP ksp, int num_eigenvalues, int argc, char* argv[])
{
    //首先用矩阵A,B创建petsc_solver
    GCGE_SOLVER *petsc_solver = GCGE_PETSC_Solver_Init(A, B, num_eigenvalues, argc, argv);
    //用户已创建好ksp,直接赋给petsc_solver
    GCGE_SOLVER_SetOpsLinearSolverWorkspace(petsc_solver, (void*)ksp);
    //将线性解法器设为KSPSolve
    petsc_solver->ops->LinearSolver = GCGE_PETSC_LinearSolver;
    return petsc_solver;
}

void GCGE_SOLVER_SetPETSCOpsLinearSolver(GCGE_SOLVER *solver, KSP ksp)
{
    GCGE_SOLVER_SetOpsLinearSolverWorkspace(solver, (void*)ksp);
    //将线性解法器设为KSPSolve
    solver->ops->LinearSolver = GCGE_PETSC_LinearSolver;
}
