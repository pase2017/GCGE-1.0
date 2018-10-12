/*
 * =====================================================================================
 *
 *       Filename:  external_petsc.c
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

#if 1
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "external_petsc.h"

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
    PetscInt           low, high, length, one = 1;
    ierr = VecGetOwnershipRange(x, &low, &high);
    length = high-low;
    ierr = VecGetArrayRead(x,&local_x);
    ierr = VecGetArrayRead(y,&local_y);
    PetscStackCallBLAS("BLASdot",*value = BLASdot_(&length,local_y,&one,local_x,&one));
    ierr = VecRestoreArrayRead(x,&local_x);
    ierr = VecRestoreArrayRead(y,&local_y);
}

#endif
