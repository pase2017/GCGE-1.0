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

/* TODO should be modified orthogonal's workspace shoud be VEC** */
#include "gcge_solver.h"

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
    const char *file_A = "../test/data/testA";
    const char *file_B = "../test/data/testB";
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
        CSR_BuildVecByMat(A, multi_vec+i);
    }
    CSR_VEC *vec;
    CSR_BuildVecByMat(A, &vec);

    for(i=0; i<vec->size; i++)
    {
       vec->Entries[i] = i+1;
    }

    int k;
    for(k=0; k<num_vec; k++)
    {
       for(i=0; i<multi_vec[k]->size; i++)
       {
	  multi_vec[k]->Entries[i] = i+k;
       }
       CSR_PrintVec(multi_vec[k]);
    }
    int start = 0, end = 2;
    /* SHOULD be modified to ops */
    GCGE_Orthogonal((void **)multi_vec, start, &end, (void *)B, 
	  solver->ops, solver->para->orth_para, solver->workspace->V_tmp, solver->workspace->subspace_dtmp);

    for(k=0; k<num_vec; k++)
    {
       CSR_PrintVec(multi_vec[k]);
    }

    double vTv;
    solver->ops->VecInnerProd((void *)multi_vec[0], (void *)multi_vec[1], &vTv);
    printf ( "vTv = %f\n", vTv );

    double mvTmv[4];
    int    mv_s[2] = {0, 0};
    int    mv_e[2] = {2, 2};
    solver->ops->MultiVecInnerProd((void **)multi_vec, (void **)multi_vec, (double *)&mvTmv, 
	  "S", mv_s, mv_e, 2, solver->ops);

    printf ( "mvTmv\n" );
    printf ( "%f\t%f\n%f\t%f\n", mvTmv[0], mvTmv[1], mvTmv[2], mvTmv[3] );
    
    double entries[6] = {0.0, 1.0, 2.0, 1.0, 2.0, 3.0};
    end = 2;
    
    //  void GCGE_OrthogonalSubspace(double *V, GCGE_INT ldV, GCGE_INT start, GCGE_INT *end, 
    //  void *B, GCGE_INT ldB, GCGE_ORTH_PARA *orth_para);
    GCGE_OrthogonalSubspace((double *)entries, 3, 0, &end,
	  NULL, -1, solver->para->orth_para);
    printf ( "entries\n" );
    printf ( "%16.5f\t%10.15f\t%10.15f\n%10.15f\t%10.15f\t%10.15f\n", 
	  entries[0], entries[1], entries[2], 
	  entries[3], entries[4], entries[5]  );


    //释放矩阵空间
    GCGE_SOLVER_Free(&solver);
    CSR_MatFree(&A);
    CSR_MatFree(&B);
    for(i=0; i<num_vec; i++)
    {
        CSR_VecFree(multi_vec+i);
    }
    free(multi_vec); multi_vec = NULL;
    CSR_VecFree(&vec);

    return 0;
}
