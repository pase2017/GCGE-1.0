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
#include "gcge_app_slepc.h"

static char help[] = "Use GCGE-SLEPc to solve an eigensystem Ax=kBx with the matrixes loaded from files.\n";

int main(int argc, char* argv[])
{
    PetscErrorCode ierr;
    SlepcInitialize(&argc,&argv,(char*)0,help);
    
    GCGE_DOUBLE t_start = GCGE_GetTime();

    //从参数中读入矩阵A, B, P的地址, 矩阵P为求解线性方程组时的预条件矩阵
    Mat A, B, P;
    char file_A[PETSC_MAX_PATH_LEN] = "fileinput";
    ierr = PetscOptionsGetString(NULL, NULL, "-mat_A", file_A, sizeof(file_A), NULL);
    char file_B[PETSC_MAX_PATH_LEN] = "fileinput";
    ierr = PetscOptionsGetString(NULL, NULL, "-mat_B", file_B, sizeof(file_B), NULL);
    char file_P[PETSC_MAX_PATH_LEN] = "fileinput";
    ierr = PetscOptionsGetString(NULL, NULL, "-mat_P", file_P, sizeof(file_P), NULL);

    //读入矩阵, 如果命令行提供了矩阵B和P的地址, 才读矩阵B和P
    SLEPC_ReadMatrixBinary(&A, file_A);
    if(strcmp(file_B, "fileinput") != 0)
    {
        SLEPC_ReadMatrixBinary(&B, file_B);
    }
    if(strcmp(file_P, "fileinput") != 0)
    {
        SLEPC_ReadMatrixBinary(&P, file_P);
    }

    //创建petsc_solver
    GCGE_SOLVER *slepc_solver;
    GCGE_SOLVER_Create(&slepc_solver);
    //设置SLEPC结构的矩阵向量操作
    GCGE_SOLVER_SetSLEPCOps(slepc_solver);

    //设置一些参数-示例
    int nev = 30;
    GCGE_SOLVER_SetNumEigen(slepc_solver, nev);//设置特征值个数
    slepc_solver->para->print_eval = 0;//设置是否打印每次迭代的特征值

    //从命令行读入GCGE_PARA中的一些参数
    GCGE_INT error = GCGE_PARA_SetFromCommandLine(slepc_solver->para, argc, argv);

    //给特征值和特征向量分配nev大小的空间
    nev = slepc_solver->para->nev;
    double *eval = (double *)calloc(nev, sizeof(double)); 
    BV evec;
    slepc_solver->ops->BuildMultiVecByMat((void*)A, (void***)(&evec), nev, slepc_solver->ops);

    //将矩阵与特征对设置到petsc_solver中
    GCGE_SOLVER_SetMatA(slepc_solver, A);
    if(strcmp(file_B, "fileinput") != 0)
    {
        GCGE_SOLVER_SetMatB(slepc_solver, B);
    }
    GCGE_SOLVER_SetEigenvalues(slepc_solver, eval);
    GCGE_SOLVER_SetEigenvectors(slepc_solver, (void**)evec);

    KSP ksp;
    //如果读入了预条件矩阵, 就使用PETSc的求解器
    if(strcmp(file_P, "fileinput") != 0)
    {
        //设定线性求解器
        SLEPC_LinearSolverCreate(&ksp, A, P);
        //PetscViewer viewer;
        //ierr = KSPView(ksp, viewer);
        //给slepc_solver设置KSP为线性求解器
        GCGE_SOLVER_SetSLEPCOpsLinearSolver(slepc_solver, ksp);
    }

    //对slepc_solver进行setup，检查参数，分配工作空间等
    GCGE_SOLVER_Setup(slepc_solver);
    //求解
    GCGE_SOLVER_Solve(slepc_solver);

    //释放特征值、特征向量、KSP空间
    free(eval); eval = NULL;
    slepc_solver->ops->FreeMultiVec((void***)(&evec), nev, slepc_solver->ops);
    if(strcmp(file_P, "fileinput") != 0)
    {
        KSPDestroy(&ksp);
    }

    //释放petsc_solver中非用户创建的空间
    GCGE_SOLVER_Free(&slepc_solver);

    //释放矩阵空间
    ierr = MatDestroy(&A);
    if(strcmp(file_B, "fileinput") != 0)
    {
        ierr = MatDestroy(&B);
    }
    if(strcmp(file_P, "fileinput") != 0)
    {
        ierr = MatDestroy(&P);
    }

    GCGE_DOUBLE t_end = GCGE_GetTime();
    GCGE_Printf("From ReadMatrix To MatDestroy, Total Time: %f\n", t_end - t_start);

    ierr = SlepcFinalize();

    return ierr;
}
