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

#include "gcge.h"
#include "gcge_app_petsc.h"
#include "external_petsc.h"

static char help[] = "Use GCGE-PETSc to solve an eigensystem Ax=kBx with the matrixes loaded from files.\n";

int main(int argc, char* argv[])
{
    PetscErrorCode ierr;
    PetscInitialize(&argc,&argv,(char*)0,help);

    //读入矩阵,其中矩阵P为求解线性方程组时的预条件矩阵
    Mat A, B;
    char file_A[PETSC_MAX_PATH_LEN] = "fileinput";
    ierr = PetscOptionsGetString(NULL, NULL, "-mat_A", file_A, sizeof(file_A), NULL);
    char file_B[PETSC_MAX_PATH_LEN] = "fileinput";
    ierr = PetscOptionsGetString(NULL, NULL, "-mat_B", file_B, sizeof(file_B), NULL);
    PETSC_ReadMatrixBinary(&A, file_A);
    PETSC_ReadMatrixBinary(&B, file_B);

    //创建petsc_solver
    GCGE_SOLVER *petsc_solver;
    GCGE_SOLVER_Create(&petsc_solver);

    //设置一些参数
    int nev = 30;
    GCGE_SOLVER_SetNumEigen(petsc_solver, nev);//设置特征值个数
    petsc_solver->para->ev_tol = 1e-12;//设置收敛准则
    petsc_solver->para->ev_max_it = 30;//设置GCGE最大迭代次数
    petsc_solver->para->print_part_time = 1;//设置是否打印每次迭代的时间统计信息

    //从命令行读入GCGE_PARA中的一些参数
    GCGE_INT error = GCGE_PARA_SetFromCommandLine(petsc_solver->para, argc, argv);

    //给特征值和特征向量分配nev大小的空间
    nev = petsc_solver->para->nev;
    double *eval = (double *)calloc(nev, sizeof(double)); 
    Vec *evec;
	Vec vector;
    ierr = MatCreateVecs(A, NULL, &vector);
    ierr = VecDuplicateVecs(vector, nev, &evec);

    //将矩阵与特征对设置到petsc_solver中
    GCGE_SOLVER_SetMatA(petsc_solver, A);
    GCGE_SOLVER_SetMatB(petsc_solver, B);
    GCGE_SOLVER_SetEigenvalues(petsc_solver, eval);
    GCGE_SOLVER_SetEigenvectors(petsc_solver, (void**)evec);

    //设置PETSC结构的矩阵向量操作
    GCGE_SOLVER_SetPETSCOps(petsc_solver);

    //设定线性求解器
    KSP ksp;
    PETSC_LinearSolverCreate(&ksp, A, A);
    PetscViewer viewer;
    ierr = KSPView(ksp, viewer);
    //给petsc_solver设置KSP为线性求解器
    GCGE_SOLVER_SetPETSCOpsLinearSolver(petsc_solver, ksp);

    //对petsc_solver进行setup，检查参数，分配工作空间等
    GCGE_SOLVER_Setup(petsc_solver);
    //求解
    GCGE_SOLVER_Solve(petsc_solver);
    //释放petsc_solver中非用户创建的空间
    GCGE_SOLVER_Free(&petsc_solver);

    //释放特征值、特征向量、KSP空间
    free(eval); eval = NULL;
    VecDestroyVecs(nev, &evec);
    KSPDestroy(&ksp);

    //释放矩阵空间
    ierr = MatDestroy(&A);
    ierr = MatDestroy(&B);

    ierr = PetscFinalize();

    return ierr;
}
