/*
 * =====================================================================================
 *
 *       Filename:  csr.c
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

#include <petscoptions.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petscviewer.h>
#include <petscksp.h>
#include <petscpc.h>

void PETSC_ReadMatrixBinary(Mat *A, const char *filename);
void PETSC_LinearSolverCreate(KSP *ksp, Mat A, Mat T);

#endif
