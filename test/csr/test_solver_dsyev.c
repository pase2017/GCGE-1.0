/*
 * =====================================================================================
 *
 *       Filename:  ex1.c
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

#include "csr.h"
#include "gcge.h"

#include "gcge_app_csr.h"

void GCGE_DenseMatEigenSolver_dsyevr(char *jobz, char *range, char *uplo, 
        GCGE_INT *nrows, GCGE_DOUBLE *a, GCGE_INT *lda, 
        GCGE_DOUBLE *vl, GCGE_DOUBLE *vu, GCGE_INT *il, GCGE_INT *iu, 
        GCGE_DOUBLE *abstol, GCGE_INT *nev, 
        GCGE_DOUBLE *eval, GCGE_DOUBLE *evec, GCGE_INT *lde, 
        GCGE_INT *isuppz, GCGE_DOUBLE *work, GCGE_INT *lwork, 
        GCGE_INT *iwork, GCGE_INT *liwork, 
        GCGE_INT *ifail, GCGE_INT *info)
{
    *lwork = 26*(*nrows);
    *liwork = 10*(*nrows);
    dsyevr_(jobz, range, uplo, nrows, a, lda, vl, vu, il, iu, 
            abstol, nev, eval, evec, lde, isuppz, work, lwork, 
            iwork, liwork, info);
}

void GCGE_SOLVER_SetOpsDsyevr(GCGE_SOLVER *solver)
{
    solver->ops->DenseMatEigenSolver = GCGE_DenseMatEigenSolver_dsyevr;
}

int main(int argc, char* argv[])
{
    GCGE_SOLVER *solver;
    int nev = 100;

    //创建矩阵
    const char *file_A = "../data/testA";
    //const char *file_B = "../data/testB";
    //const char *file_A = "../data/A_5.txt";
    //const char *file_B = "../data/M_5.txt";
    //const char *file_A = "/home/nzhang/A_7.txt";
    //const char *file_B = "/home/nzhang/M_7.txt";
//    const char *file_A = "/home/nzhang/Andrews.csr";
    //const char *file_B = "/home/nzhang/M_7.txt";
    CSR_MAT *A = CSR_ReadMatFile(file_A);
    //CSR_MAT *B = CSR_ReadMatFile(file_B);

    //创建特征对空间
    double *eval = (double *)calloc(nev, sizeof(double));
    CSR_VEC **evec = (CSR_VEC**)malloc(nev*sizeof(CSR_VEC*));
    int i = 0; 
    for(i=0; i<nev; i++)
    {
        CSR_BuildVecByMat(A, evec+i);
    }

    //------------------- dsyevx -----------------------
    double tx1 = GCGE_GetTime();
    srand(1);
    GCGE_INT help = GCGE_SOLVER_Create(&solver, argc, argv);
    if(help)
    {
        GCGE_SOLVER_Free(&solver);
        exit(0);
    }
    GCGE_SOLVER_SetEigenvalues(solver, eval);
    GCGE_SOLVER_SetEigenvectors(solver, (void**)evec);
    GCGE_SOLVER_SetMatA(solver, A);
    //GCGE_SOLVER_SetMatB(solver, B);
    GCGE_SOLVER_SetCSROps(solver);
    GCGE_SOLVER_SetNumEigen(solver, nev);
    //solver->para->ev_max_it = 3;
    solver->para->ev_tol = 1e-4;
    //solver->para->orth_para->print_orth_zero= 1;
    //solver->para->orth_para->orth_zero_tol = 1e-3;
    GCGE_SOLVER_Setup(solver);
    GCGE_SOLVER_Solve(solver);
    GCGE_SOLVER_Free(&solver);
    double tx2 = GCGE_GetTime();

    //------------------- dsyevr -----------------------
    //srand((unsigned)time(NULL));
    srand(1);
    double tr1 = GCGE_GetTime();
    help = GCGE_SOLVER_Create(&solver, argc, argv);
    if(help)
    {
        GCGE_SOLVER_Free(&solver);
        exit(0);
    }
    GCGE_SOLVER_SetEigenvalues(solver, eval);
    GCGE_SOLVER_SetEigenvectors(solver, (void**)evec);
    GCGE_SOLVER_SetMatA(solver, A);
    //GCGE_SOLVER_SetMatB(solver, B);
    GCGE_SOLVER_SetCSROps(solver);
    GCGE_SOLVER_SetOpsDsyevr(solver);
    GCGE_SOLVER_SetNumEigen(solver, nev);
    solver->para->ev_tol = 1e-4;
    GCGE_SOLVER_Setup(solver);
    GCGE_SOLVER_Solve(solver);
    GCGE_SOLVER_Free(&solver);
    double tr2 = GCGE_GetTime();
    printf("dsyevx time: %f, dsyevr time: %f\n", tx2-tx1, tr2-tr1);

    //释放矩阵空间
    CSR_MatFree(&A);
    //CSR_MatFree(&B);
    for(i=0; i<nev; i++)
    {
        CSR_VecFree(evec+i);
    }
    free(evec); evec = NULL;
    free(eval); eval = NULL;
    //释放特征对空间

    return 0;
}

