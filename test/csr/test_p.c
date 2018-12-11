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
    CSR_VEC **multi_vec;
    solver->ops->MultiVecCreateByMat((void***)(&multi_vec), num_vec, (void*)A, solver->ops);

    int num_vec_2 = 2;
    CSR_VEC **multi_vec_2;
    solver->ops->MultiVecCreateByMat((void***)(&multi_vec_2), num_vec_2, (void*)A, solver->ops);

    int j = 0;
    solver->ops->MultiVecSetRandomValue((void**)multi_vec, num_vec, solver->ops);

    int k;
    printf("before P:\n");
    solver->ops->MultiVecPrint((void**)multi_vec, num_vec);

    int dim_x = 2;
    solver->workspace->dim_xpw = num_vec;
    solver->workspace->dim_x = dim_x;
    solver->workspace->unconv_bs = 1;
    solver->workspace->last_dim_x = 1;
    solver->workspace->unlock[0] = 0;
    //solver->para->cg_max_it = 3;
    //solver->para->cg_rate = 1e-10;

    double subspace_evec[9] = {5.0, 12.0, 3.0, 34.0, 5.0, 6.0, 1.0, 1.0, 1.0};

    //对基底multi_vec正交化
    GCGE_Orthonormalization((void**)multi_vec, 0, &num_vec, NULL, solver->ops, solver->para->orth_para, 
            solver->workspace->V_tmp, solver->workspace->subspace_dtmp);

    printf("after  X orth:\n");
    solver->ops->MultiVecPrint((void**)multi_vec, num_vec);

    GCGE_OrthonormalizationSubspace(subspace_evec, num_vec, 0, &dim_x, NULL, -1, solver->para->orth_para);
    //基底multi_vec线性组合得到multi_vec_2, dim_x个
    int mv_s[2] = { 0, 0 };
    int mv_e[2] = { num_vec, dim_x };
    solver->ops->MultiVecLinearComb((void**)multi_vec, (void**)multi_vec_2,
            (int*)mv_s, (int*)mv_e, subspace_evec, num_vec, NULL, -1, solver->ops);
    /*
    printf("b, subspace_evec: \n");
    printf("%f\t%f\t%f\n", subspace_evec[0], subspace_evec[3], subspace_evec[6]);
    printf("%f\t%f\t%f\n", subspace_evec[1], subspace_evec[4], subspace_evec[7]);
    */

    GCGE_ComputeP((double*)subspace_evec, (void **)multi_vec,
        solver->ops, solver->para, solver->workspace);
    printf("subspace_evec: \n");
    printf("%f\t%f\t%f\n", subspace_evec[0], subspace_evec[3], subspace_evec[6]);
    printf("%f\t%f\t%f\n", subspace_evec[1], subspace_evec[4], subspace_evec[7]);
    printf("%f\t%f\t%f\n", subspace_evec[2], subspace_evec[5], subspace_evec[8]);

    mv_s[0] = 0;
    mv_e[0] = dim_x;
    mv_s[1] = 0;
    mv_e[1] = dim_x;
    solver->ops->MultiVecAxpby(1.0, (void**)multi_vec_2, 0.0, (void**)multi_vec, mv_s, mv_e, solver->ops);
    int i = 0;

    printf("after  P orth:\n");
    for(k=0; k<num_vec; k++)
    {
       CSR_PrintVec(multi_vec[k]);
    }

    double alpha[9];
    mv_s[0] = 0;
    mv_e[0] = 3;
    mv_s[1] = 0;
    mv_e[1] = 3;
    solver->ops->MultiVecInnerProd((void**)multi_vec, (void**)multi_vec, 
            (double*)alpha, "nonsym", mv_s, mv_e, num_vec, solver->ops);
    printf("x^Tp: \n");
    printf("%f\t%f\t%f\n", alpha[0], alpha[3], alpha[6]);
    printf("%f\t%f\t%f\n", alpha[1], alpha[4], alpha[7]);
    printf("%f\t%f\t%f\n", alpha[2], alpha[5], alpha[8]);
    
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
        CSR_VecDestroy(multi_vec+i);
    }
    free(multi_vec); multi_vec = NULL;
    for(i=0; i<num_vec_2; i++)
    {
        CSR_VecDestroy(multi_vec_2+i);
    }
    free(multi_vec_2); multi_vec_2 = NULL;

    return 0;
}
