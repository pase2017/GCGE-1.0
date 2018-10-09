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
    //创建矩阵
    Mat A, B;
    PetscErrorCode ierr;

    PetscInitialize(&argc,&argv,(char*)0,help);
    double t1 = 0.0;
    double t2 = 0.0;

    const char *file_A = "../data/A_5.petsc.bin";
    const char *file_B = "../data/M_5.petsc.bin";
    PETSC_ReadMatrixBinary(&A, file_A);
    PETSC_ReadMatrixBinary(&B, file_B);
    int nev = 30;
    GCGE_SOLVER *petsc_solver = GCGE_PETSC_Solver_Init_KSPLinearSolver(A, B, A, nev, argc,  argv);   

    petsc_solver->para->ev_tol = 1e-12;
    petsc_solver->para->orth_para->print_orth_zero = 1;
    petsc_solver->para->ev_max_it = 30;
    petsc_solver->para->print_part_time = 1;

    GCGE_SOLVER_Solve(petsc_solver);  

    GCGE_SOLVER_Free_All(&petsc_solver);
    //释放矩阵空间
    ierr = MatDestroy(&A);
    ierr = MatDestroy(&B);
    ierr = PetscFinalize();

    return ierr;
}
