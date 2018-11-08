/*
 * =====================================================================================
 *
 *       Filename:  test_orthogonal.c
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

    int num_vec = 3;
    CSR_VEC **multi_vec = (CSR_VEC**)malloc(num_vec*sizeof(CSR_VEC*));
    int i = 0; 
    for(i=0; i<num_vec; i++)
    {
        CSR_BuildVecByMat(A, multi_vec+i);
    }

    int num_vec_2 = 2;
    CSR_VEC **multi_vec_2 = (CSR_VEC**)malloc(num_vec*sizeof(CSR_VEC*));
    for(i=0; i<num_vec_2; i++)
    {
        CSR_BuildVecByMat(A, multi_vec_2+i);
    }

    int j = 0;
    solver->ops->MultiVecSetRandomValue((void**)multi_vec, num_vec, solver->ops);

    int k;
    printf("before P:\n");
    for(k=0; k<num_vec; k++)
    {
       CSR_PrintVec(multi_vec[k]);
    }
    double a[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    solver->workspace->dim_x = 1;
    solver->workspace->dim_xpw = 3;
    solver->workspace->last_dim_x = 1;
    solver->workspace->unconv_bs = 1;
    solver->workspace->unlock[0] = 0;
    //solver->para->cg_max_it = 3;
    //solver->para->cg_rate = 1e-10;

    double eval[2] = {1.0, 2.0};
    double subspace_evec[9] = {1.0/sqrt(14.0), 2.0/sqrt(14.0), 3.0/sqrt(14.0), 4.0, 5.0, 6.0, 1.0, 1.0, 1.0};

    int start[2] = {0, 0};
    int end[2] = {3, 1};
    int lda = 3;
    GCGE_Orthogonal((void**)multi_vec, 0, &lda, NULL, solver->ops, solver->para->orth_para, solver->workspace);
    solver->ops->MultiVecLinearComb((void**)multi_vec, (void**)multi_vec_2,
            (int*)start, (int*)end, subspace_evec, lda, NULL, -1, solver->ops);

    GCGE_ComputeP((double*)subspace_evec, (void **)multi_vec,
        solver->ops, solver->para, solver->workspace);

    printf("after  P:\n");
    for(k=0; k<num_vec; k++)
    {
       CSR_PrintVec(multi_vec[k]);
    }

    double alpha;
    solver->ops->VecInnerProd(multi_vec_2[0], multi_vec[1], &alpha);
    printf("x^Tp: %f\n", alpha);
    
    /*
    for(i=0; i<2; i++)
    {
        solver->ops->MatDotVec((void*)A, multi_vec[2+i], multi_vec_2[0]);
        solver->ops->MatDotVec((void*)B, multi_vec[i],   multi_vec_2[1]);
        solver->ops->VecAxpby(1.0, multi_vec_2[0], -eval[i], multi_vec_2[1]);
        solver->ops->VecInnerProd(multi_vec_2[1], multi_vec_2[1], &alpha);
        printf("A y - lambda B x: norm[%d]: %f\n", i, sqrt(alpha));
    }
    */

    //释放矩阵空间
    GCGE_SOLVER_Free(&solver);
    CSR_MatFree(&A);
    CSR_MatFree(&B);
    for(i=0; i<num_vec; i++)
    {
        CSR_VecFree(multi_vec+i);
    }
    free(multi_vec); multi_vec = NULL;
    for(i=0; i<num_vec_2; i++)
    {
        CSR_VecFree(multi_vec_2+i);
    }
    free(multi_vec_2); multi_vec_2 = NULL;

    return 0;
}
