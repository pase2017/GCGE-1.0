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

// W: A*W = lambda * B * X
// res = A*W-lambda*B*X
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
    CSR_VEC **multi_vec;
    solver->ops->BuildMultiVecByMat((void*)A, (void***)(&multi_vec), num_vec, solver->ops);
    solver->ops->MultiVecSetRandomValue((void**)multi_vec, num_vec, solver->ops);
    solver->ops->PrintMultiVec((void**)multi_vec, 3);

    int num_vec_2 = 3;
    CSR_VEC **multi_vec_2;
    solver->ops->BuildMultiVecByMat((void*)A, (void***)(&multi_vec_2), num_vec_2, solver->ops);

    solver->workspace->dim_xp = 1;
    solver->workspace->unconv_bs = 1;
    solver->workspace->unlock[0] = 0;

    double eval[3] = {1.0, 2.0, 3.0};
    GCGE_ComputeW((void*)A, (void*)B, (void**)multi_vec, eval,
            solver->ops, solver->para, solver->workspace);

    int i = 0;
    //W=multi_vec[1]
    //X=multi_vec[0]
    //AW=multi_vec_2[1]
    //BX=multi_vec_2[0]
    //multi_vec_2[0]= AW-lambda*BX
    solver->ops->MatDotVec((void*)A, (void*)(multi_vec[1]), (void*)(multi_vec_2[1]));
    solver->ops->MatDotVec((void*)B, (void*)(multi_vec[0]), (void*)(multi_vec_2[0]));
    printf("Aw,Bx");
    solver->ops->PrintMultiVec((void**)multi_vec_2, 2);
    solver->ops->VecAxpby(1.0, (void*)multi_vec_2[1], -eval[0], (void*)multi_vec_2[0]);

    printf("Aw-lambdaBx");
    solver->ops->PrintMultiVec((void**)multi_vec_2, 2);

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
