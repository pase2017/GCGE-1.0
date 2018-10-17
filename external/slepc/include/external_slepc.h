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

#include <slepcbv.h>

void SLEPC_ReadMatrixBinary(Mat *A, const char *filename);
void SLEPC_LinearSolverCreate(KSP *ksp, Mat A, Mat T);
void SLEPC_VecLocalInnerProd(Vec x, Vec y, double *value);
#endif
