/*
 * =====================================================================================
 *
 *       Filename:  test_solver.c
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

//#include "memwatch.h"
#include <petscsys.h>
#include <petscviewer.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscpc.h>

#include "gcge.h"
#include "gcge_app_petsc.h"
#include "external_petsc.h"

static char help[] = "Use GCGE-PETSc to solve an eigensystem Ax=kBx with the matrixes loaded from files.\n";

#if 0
void LinearSolverCreate(KSP *ksp, Mat A, Mat T)
{
    PetscErrorCode ierr;
    ierr = KSPCreate(PETSC_COMM_WORLD,ksp);
    ierr = KSPSetOperators(*ksp,A,T);
    ierr = KSPSetFromOptions(*ksp);
    ierr = KSPSetInitialGuessNonzero(*ksp, 1);
}

void GCGE_OPS_SetLinearSolverWorkspace(GCGE_OPS *ops, void *linear_solver_workspace)
{
    ops->linear_solver_workspace = linear_solver_workspace;
}

void GCGE_SOLVER_SetOpsLinearSolverWorkspace(GCGE_SOLVER *solver, void *linear_solver_workspace)
{
    GCGE_OPS_SetLinearSolverWorkspace(solver->ops, linear_solver_workspace);
}

void GCGE_SOLVER_SetKSPLinearSolver(GCGE_SOLVER *solver, void *A, void *P)
{
    KSP ksp;
    LinearSolverCreate(&ksp, (Mat)A, (Mat)P);
    GCGE_SOLVER_SetOpsLinearSolverWorkspace(solver, (void*)ksp);
    ops->LinearSolver = GCGE_PETSC_LinearSolver;
}
#endif
void GCGE_PETSC_KSPLinearSolver(void *Matrix, void *b, void *x, struct GCGE_OPS_ *ops)
{
    PetscErrorCode ierr;
    ierr = KSPSolve((KSP)(ops->linear_solver_workspace), (Vec)b, (Vec)x);
}


int main(int argc, char* argv[])
{
    //创建矩阵
    Mat A;
    PetscErrorCode ierr;

    PetscInitialize(&argc,&argv,(char*)0,help);
    double t1 = 0.0;
    double t2 = 0.0;

    const char *file_A = "../data/A_5.petsc.bin";
    PETSC_ReadMatrixBinary(&A, file_A);

    GCGE_OPS *ops;
    GCGE_OPS_Create(&ops);
    GCGE_PETSC_SetOps(ops);
    GCGE_INT error = GCGE_OPS_Setup(ops);

    KSP ksp;
    PETSC_LinearSolverCreate(&ksp, A, A);

    //添加下面这几行设置后可以使用KSP精确求解线性方程组
#if 0
    ierr = KSPSetType(ksp, KSPCG);
    ierr = KSPSetTolerances(ksp, 1e-15, PETSC_DEFAULT, PETSC_DEFAULT, 1000);
    PC pc;
    ierr = KSPGetPC(ksp, &pc);
    ierr = PCSetType(pc, PCHYPRE);
    ierr = PCHYPRESetType(pc, "boomeramg");
#endif

    GCGE_OPS_SetLinearSolverWorkspace(ops, (void*)ksp);
    ops->LinearSolver = GCGE_PETSC_KSPLinearSolver;

    Vec x;
    Vec rhs;
    Vec tmp;
    ops->BuildVecByMat((void*)A, (void*)(&x));
    ops->BuildVecByMat((void*)A, (void*)(&rhs));
    ops->BuildVecByMat((void*)A, (void*)(&tmp));
    ops->VecSetRandomValue((void*)x);
    ops->VecSetRandomValue((void*)rhs);

    ops->LinearSolver((void*)A, (void*)rhs, (void*)x, ops);

    ops->MatDotVec((void*)A, (void*)x, (void*)tmp);
    ops->VecAxpby(1.0, (void*)rhs, -1.0, (void*)tmp);
    double norm = GCGE_VecNorm((void*)tmp, ops);

    GCGE_Printf("|| A * x - rhs ||_2: %e\n", norm);

    ops->FreeVec((void*)(&x));
    ops->FreeVec((void*)(&rhs));
    ops->FreeVec((void*)(&tmp));

    GCGE_OPS_Free(&ops);
    //释放矩阵空间
    ierr = MatDestroy(&A);

    ierr = PetscFinalize();

    return ierr;
}
