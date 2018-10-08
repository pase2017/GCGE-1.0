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

PetscErrorCode ReadPetscMatrixBinary(Mat *A, const char *filename)
{
    PetscErrorCode ierr;
    PetscViewer    viewer;
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Getting matrix...\n"); CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer); CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, A); CHKERRQ(ierr);
    ierr = MatSetFromOptions(*A); CHKERRQ(ierr);
    ierr = MatLoad(*A, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Getting matrix... Done\n"); CHKERRQ(ierr);
    return ierr;
}

#endif
