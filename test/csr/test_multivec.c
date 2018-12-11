/*
 * =====================================================================================
 *
 *       Filename:  test_ops.c
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

//#include "memwatch.h"

#include "gcge_app_csr.h"

int main(int argc, char* argv[])
{
//    mwInit();
    srand((unsigned)time(NULL));

    //创建矩阵
    const char *file_A = "../data/testA";
    CSR_MAT *A = CSR_ReadMatFile(file_A);

    GCGE_OPS *ops = NULL;
    GCGE_OPS_Create(&ops);
    GCGE_CSR_SetOps(ops);
    GCGE_OPS_Setup(ops);

    int num_vec_x = 2;
    int num_vec_y = 3;
    CSR_VEC **multi_vec_x;
    CSR_VEC **multi_vec_y;
    printf ( "BuildMultiVecByMat\n" );
    ops->BuildMultiVecByMat((void *)A, (void ***)(&multi_vec_x), num_vec_x, ops);
    ops->BuildMultiVecByMat((void *)A, (void ***)(&multi_vec_y), num_vec_y, ops);

    int i, k;
    for(k=0; k<num_vec_x; k++)
    {
       for(i=0; i<multi_vec_x[k]->size; i++)
       {
	  multi_vec_x[k]->Entries[i] = i+k;
       }
    }
    for(k=0; k<num_vec_y; k++)
    {
       for(i=0; i<multi_vec_y[k]->size; i++)
       {
	  multi_vec_y[k]->Entries[i] = i+k-1;
       }
    }

    printf ( "x\n" );
    ops->PrintMultiVec((void **)multi_vec_x, num_vec_x);
    printf ( "y\n" );
    ops->PrintMultiVec((void **)multi_vec_y, num_vec_y);


    double mvTmv[4], norm[2];
    int    mv_s[2] = {0, 1};
    int    mv_e[2] = {2, 3};
//    printf ( "MultiVecSetRandomValue x, y\n" );
//    ops->MultiVecSetRandomValue((void **)multi_vec_x, mv_s, mv_e, ops);
//    ops->MultiVecSetRandomValue((void **)multi_vec_y, mv_s+1, mv_e+1, ops);
//    printf ( "x\n" );
//    ops->PrintMultiVec((void **)multi_vec_x, num_vec_x);
//    printf ( "y\n" );
//    ops->PrintMultiVec((void **)multi_vec_y, num_vec_y);

    printf ( "MultiVecNorm \n" );
    printf ( "x norm\n" );
    ops->MultiVecNorm((void **)multi_vec_x, norm, mv_s, mv_e, ops);
    printf ( "%f\t%f\n", norm[0], norm[1] );
    printf ( "y norm\n" );
    ops->MultiVecNorm((void **)multi_vec_y, norm, mv_s+1, mv_e+1, ops);
    printf ( "%f\t%f\n", norm[0], norm[1] );

    printf ( "MatDotMultiVec y = Ax\n" );
    printf ( "A\n" );
    CSR_PrintMat(A);
    ops->MatDotMultiVec((void *)A, (void **)multi_vec_x, (void **)multi_vec_y, 
	  mv_s, mv_e, ops);
    printf ( "y\n" );
    ops->PrintMultiVec((void **)multi_vec_y, num_vec_y);

    double alpha = 2.0, beta = -1.0;
    printf ( "MultiVecAxpby y = %0.4fx+%0.4fy\n", alpha, beta );
    ops->MultiVecAxpby(alpha, (void **)multi_vec_x, beta, (void **)multi_vec_y, 
	  mv_s, mv_e, ops);
    ops->PrintMultiVec((void **)multi_vec_y, num_vec_y);

    printf ( "MultiVecLinearComb y = x coord\n" );
    int ldcoord = 2;
    double coord[6] = {1.0, 1.0, -1.0, -1.0, 1.0, -1.0};
    printf ( "x\n" );
    ops->PrintMultiVec((void **)multi_vec_x, num_vec_x);
    printf ( "coord\n" );
    printf ( "%f\t%f\t%f\n", coord[0], coord[2], coord[4] );
    printf ( "%f\t%f\t%f\n", coord[1], coord[3], coord[5] );

    mv_s[0] = 0; mv_s[1] = 0;
    mv_e[0] = 2; mv_e[1] = 3;
    ops->MultiVecLinearComb((void **)multi_vec_x, (void **)multi_vec_y, 
	  mv_s, mv_e, coord, ldcoord, NULL, -1, ops);
    printf ( "y\n" );
    ops->PrintMultiVec((void **)multi_vec_y, num_vec_y);

    printf ( "MultiVecInnerProd xTy\n" );
    double xTy[6] = {1.0, 1.0, -1.0, -1.0, 1.0, -1.0};
    ops->MultiVecInnerProd((void **)multi_vec_x, (void **)multi_vec_y, 
	  (double *)&xTy, "n", mv_s, mv_e, num_vec_x, ops);
    printf ( "xTy\n" );
    printf ( "%f\t%f\t%f\n", xTy[0], xTy[2], xTy[4] );
    printf ( "%f\t%f\t%f\n", xTy[1], xTy[3], xTy[5] );

    printf ( "MultiVecSwap x<->y\n" );
    mv_s[0] = 0; mv_s[1] = 1;
    mv_e[0] = 2; mv_e[1] = 3;
    ops->MultiVecSwap((void **)multi_vec_x, (void **)multi_vec_y, 
	  mv_s, mv_e, ops);
    printf ( "x\n" );
    ops->PrintMultiVec((void **)multi_vec_x, num_vec_x);
    printf ( "y\n" );
    ops->PrintMultiVec((void **)multi_vec_y, num_vec_y);
    printf ( "Get and Restor x\n" );
    CSR_VEC *tmp;
    ops->GetVecFromMultiVec((void **)multi_vec_x, 1, (void *)&tmp);
    tmp->Entries[0] = 100.0;
    ops->RestoreVecForMultiVec((void **)multi_vec_x, 1, (void *)&tmp);
    printf ( "x\n" );
    ops->PrintMultiVec((void **)multi_vec_x, num_vec_x);


    //释放矩阵空间
    CSR_MatFree(&A);

    ops->MultiVecDestroy ((void ***)&multi_vec_x, num_vec_x, ops);
    ops->MultiVecDestroy ((void ***)&multi_vec_y, num_vec_y, ops);

    GCGE_OPS_Free(&ops);
//    mwTerm();
    return 0;
}
