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
    const char *file_A = "../data/A_3.txt";
    const char *file_B = "../data/M_3.txt";
    CSR_MAT *A = CSR_ReadMatFile(file_A);
    CSR_MAT *B = CSR_ReadMatFile(file_B);

//    CSR_PrintMat(A);
//    CSR_PrintMat(B);

    solver->para->nev = 2;
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
        CSR_VecCreateByMat(multi_vec+i, A);
    }
//    CSR_VEC *vec;
//    CSR_VecCreateByMat(A, &vec);
//
//    for(i=0; i<vec->size; i++)
//    {
//       vec->Entries[i] = i+1;
//    }

    int k;
    /*
    printf ( "A\n" );
    CSR_PrintMat(A);
    printf ( "V\n" );
    for(k=0; k<num_vec; k++)
    {
       for(i=0; i<multi_vec[k]->size; i++)
       {
	  multi_vec[k]->Entries[i] = i+k+3;
       }
       CSR_PrintVec(multi_vec[k]);
    }
    */

    CSR_VEC **multi_vec_tmp;
    solver->ops->MultiVecCreateByMat((void ***)&multi_vec_tmp, 3, (void *)A, solver->ops);

    double vTAw[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

    GCGE_ComputeSubspaceMatrixVTAW((void **)multi_vec, (void *)A, 0, 3,  
	  (double *)vTAw, solver->ops, (void **)multi_vec_tmp);
    printf ( "VTAW\n" );
    printf ( "%f\t%f\t%f\n", vTAw[0], vTAw[3], vTAw[6] );
    printf ( "%f\t%f\t%f\n", vTAw[1], vTAw[4], vTAw[7] );
    printf ( "%f\t%f\t%f\n", vTAw[2], vTAw[5], vTAw[8] );

    double *p = solver->workspace->subspace_matrix;
    for (i = 0; i < 9; ++i)
    {
       p[i] = (double)rand()/((double)RAND_MAX+1);
       p[i] = 0;
    }
    p[0] = 1;
    p[4] = 2;
    p[8] = 3;

    printf ( "%f\t%f\t%f\n", p[0], p[3], p[6] );
    printf ( "%f\t%f\t%f\n", p[1], p[4], p[7] );
    printf ( "%f\t%f\t%f\n", p[2], p[5], p[8] );
 
    solver->workspace->dim_xp  = 1;
    solver->workspace->dim_xpw = 3;
    solver->workspace->last_dim_xpw = 3;

    double *q = solver->workspace->subspace_evec;
    for (i = 0; i < 3; ++i)
    {
       q[i] = (double)rand()/((double)RAND_MAX+1);
       printf ( "%f\n", q[i] );
    }
 

//    GCGE_ComputeSubspaceMatrix((void *)A, (void **)multi_vec, 
//	  solver->ops,  solver->para, solver->workspace);

    
//    CSR_PrintVec(solver->workspace->V_tmp[8]);
    printf ( "%f\t%f\t%f\n", p[0], p[3], p[6] );
    printf ( "%f\t%f\t%f\n", p[1], p[4], p[7] );
    printf ( "%f\t%f\t%f\n", p[2], p[5], p[8] );

//    GCGE_ComputeSubspaceEigenpairs(solver->workspace->subspace_matrix,  
//	  solver->workspace->eval, solver->workspace->subspace_evec,  
//	  solver->ops, solver->para, solver->workspace);
//    printf ( "%f, %f, %f\n",  solver->workspace->eval[0], 
//	   solver->workspace->eval[1],
//	    solver->workspace->eval[2] );


    double eval[3] = {-1, -1, -1};
    printf ( "--------------------------------\n" );
    //CSR_PrintMat(A);
    //CSR_PrintMat(B);
    GCGE_Solve((void *)A, (void *)B, (GCGE_DOUBLE *)eval, (void **)multi_vec, 
	  solver->para, solver->ops, solver->workspace);

    printf ( "%f, %f, %f\n",  eval[0], eval[1], eval[2] );
    printf ( "--------------------------------\n" );
    //CSR_PrintMat(A);

    //释放矩阵空间
    solver->ops->MultiVecDestroy((void ***)&multi_vec_tmp, 3, solver->ops);

    GCGE_SOLVER_Free(&solver);
    CSR_MatFree(&A);
    CSR_MatFree(&B);
    for(i=0; i<num_vec; i++)
    {
        CSR_VecDestroy(multi_vec+i);
    }
    free(multi_vec); multi_vec = NULL;
//    CSR_VecDestroy(&vec);

    return 0;
}
