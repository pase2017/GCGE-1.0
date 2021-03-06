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
#include "csr.h"
#include "gcge.h"
#include "gcge_app_csr.h"

int main(int argc, char* argv[])
{
    //    mwInit();
    srand((unsigned)time(NULL));
    //创建
    GCGE_SOLVER *solver;
    GCGE_INT help = GCGE_SOLVER_Create(&solver, argc, argv);
    if(help)
    {
        GCGE_SOLVER_Free(&solver);
        return 0;
    }

    //创建矩阵
    //const char *file_A = "../data/testA";
    //const char *file_B = "../data/testB";
    const char *file_A = "../data/A_5.txt";
    const char *file_B = "../data/M_5.txt";
    CSR_MAT *A = CSR_ReadMatFile(file_A);
    CSR_MAT *B = CSR_ReadMatFile(file_B);

    /* TODO should be modidfied */
    //    int nev = GCGE_SOLVER_GetNumEigen(solver);
    //    GCGE_SOLVER_SetNumEigen(solver, nev);
    solver->para->nev = 10;
    solver->para->ev_tol = 1e-8;
    solver->para->orth_para->print_orth_zero = 1;
    solver->para->dirichlet_boundary = 1;
    int nev = solver->para->nev;
    double *eval = (double *)calloc(nev, sizeof(double));
    CSR_VEC **evec = (CSR_VEC**)malloc(nev*sizeof(CSR_VEC*));
    solver->para->cg_max_it = 10;
    solver->para->ev_max_it = 30;
    int i = 0; 
    //这里为什么不一次性生成 CSR_MultiVecCreateByMat?
    for(i=0; i<nev; i++)
    {
        CSR_VecCreateByMat(evec+i, A);
    }

    GCGE_SOLVER_SetMatA(solver, A);
    GCGE_SOLVER_SetMatB(solver, B);
    GCGE_SOLVER_SetCSROps(solver);
    GCGE_SOLVER_SetEigenvalues(solver, eval);
    GCGE_SOLVER_SetEigenvectors(solver, (void**)evec);
    /* TODO  */
    //    GCGE_SOLVER_SetEigenpaires(solver, eval);

    //setup and solve
    GCGE_SOLVER_Setup(solver);
    printf("line 80\n");
    GCGE_SOLVER_Solve(solver);
    GCGE_SOLVER_Free(&solver);

    //printf ( "A\n" );
    //CSR_PrintMat(A);
    //printf ( "B\n" );
    //CSR_PrintMat(B);
    //printf ( "lambda\n" );
    //for(i=0; i<nev; i++)
    //{
    //    printf("eval[%d] = %18.15lf\n", i, eval[i]);
    //}
    //printf ( "evec\n" );
    //for(i=0; i<nev; i++)
    //{
    //    CSR_PrintVec(evec[i]);
    //}


    //释放矩阵空间
    CSR_MatFree(&A);
    CSR_MatFree(&B);
    for(i=0; i<nev; i++)
    {
        CSR_VecDestroy(evec+i);
    }
    free(evec); evec = NULL;
    free(eval); eval = NULL;
    //释放特征对空间

    //    mwTerm();
    return 0;
}
