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

static char help[] = "Use GCGE-PETSc to solve an eigensystem Ax=kBx with the matrixes loaded from files.\nKSP is given by user\n";

int main(int argc, char* argv[])
{
    //创建矩阵
    Mat A, B, P;
    PetscErrorCode ierr;

    PetscInitialize(&argc,&argv,(char*)0,help);

    //读入矩阵
    const char *file_A = "../data/A_5.petsc.bin";
    const char *file_B = "../data/M_5.petsc.bin";
    const char *file_P = "../data/A_5.petsc.bin";
    PETSC_ReadMatrixBinary(&A, file_A);
    PETSC_ReadMatrixBinary(&B, file_B);
    PETSC_ReadMatrixBinary(&P, file_P);
    //如果用户需要从命令行读入特征值个数，这里nev需要设置为-1
    int nev = 30;

    //设定线性求解器
    KSP ksp;
    PETSC_LinearSolverCreate(&ksp, A, P);
    PetscViewer viewer;
    ierr = KSPView(ksp, viewer);

    GCGE_SOLVER *petsc_solver;
    petsc_solver = GCGE_PETSC_Solver_Init_KSPGivenByUser(A, B, ksp, nev, argc, argv);

    petsc_solver->para->ev_tol = 1e-12;
    petsc_solver->para->orth_para->print_orth_zero = 1;
    petsc_solver->para->ev_max_it = 30;
    petsc_solver->para->block_size = 15;
    petsc_solver->para->print_part_time = 1;

    GCGE_SOLVER_Solve(petsc_solver);  
    double *eval;
    Vec    *evec;
    GCGE_SOLVER_GetEigenvalues(petsc_solver, &eval);
    GCGE_SOLVER_GetEigenvectors(petsc_solver,(void***)(&evec));
    printf("eval[0]: %e\n", eval[0]);

    GCGE_SOLVER_Free_All(&petsc_solver);
    //释放矩阵空间
    ierr = MatDestroy(&A);
    ierr = MatDestroy(&B);
    ierr = PetscFinalize();

    return ierr;
}
