/*
 * =====================================================================================
 *
 *       Filename:  test_matvec.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/27/2018 04:47:29 PM
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

int main(int argc, char* argv[])
{
    srand((unsigned)time(NULL));
    GCGE_SOLVER *solver;
    GCGE_INT error = GCGE_SOLVER_Create(&solver, argc, argv);
    if(error)
    {
        GCGE_SOLVER_Free(&solver);
        exit(0);
    }

    //创建矩阵
    const char *file_A = "../data/testA";
    const char *file_B = "../data/testB";
    CSR_MAT *A = CSR_ReadMatFile(file_A);
    CSR_MAT *B = CSR_ReadMatFile(file_B);

    CSR_PrintMat(A);
    CSR_PrintMat(B);

    //setup
    GCGE_SOLVER_SetMatA(solver, A);
    GCGE_SOLVER_SetMatB(solver, B);
    GCGE_SOLVER_SetCSROps(solver);
    GCGE_SOLVER_Setup(solver);

    int num_vec = 2;
    CSR_VEC **multi_vec = (CSR_VEC**)malloc(num_vec*sizeof(CSR_VEC*));
    int i = 0; 
    for(i=0; i<num_vec; i++)
    {
        CSR_VecCreateByMat(A, multi_vec+i);
    }
    CSR_VEC *vec;
    CSR_VecCreateByMat(A, &vec);

    for(i=0; i<vec->size; i++)
    {
       vec->Entries[i] = i+1;
    }

    CSR_VEC *tmp;
    CSR_VecCreateByMat(A, &tmp);

    double vAv = GCGE_VecMatrixVec((void *)vec, (void *)A, (void *)vec, (void *)tmp, solver->ops);
    double vn = GCGE_VecNorm((void *)vec, solver->ops);

    printf ( "vAv = %f\n", vAv );
    printf ( "vn = %f\n", vn );



    //释放矩阵空间
    GCGE_SOLVER_Free(&solver);
    CSR_MatFree(&A);
    CSR_MatFree(&B);
    for(i=0; i<num_vec; i++)
    {
        CSR_VecDestroy(multi_vec+i);
    }
    free(multi_vec); multi_vec = NULL;
    CSR_VecDestroy(&vec);
    CSR_VecDestroy(&tmp);

    return 0;
}
