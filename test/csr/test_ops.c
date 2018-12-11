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
    const char *file_B = "../data/testB";
    CSR_MAT *A = CSR_ReadMatFile(file_A);
    CSR_MAT *B = CSR_ReadMatFile(file_B);

    GCGE_OPS *ops;
    GCGE_OPS_Create(&ops);
    GCGE_CSR_SetOps(ops);
    GCGE_OPS_Setup(ops);

    int num_vec_x = 2;
    int num_vec_y = 3;
    CSR_VEC **multi_vec_x;
    CSR_VEC **multi_vec_y;
    printf ( "MultiVecCreateByMat\n" );
    ops->MultiVecCreateByMat((void ***)(&multi_vec_x), num_vec_x, (void *)A, ops);
    ops->MultiVecCreateByMat((void ***)(&multi_vec_y), num_vec_y, (void *)A, ops);

//    printf ( "MultiVecSetRandomValue x, y\n" );
//    ops->MultiVecSetRandomValue((void **)multi_vec_x, num_vec_x, ops);
//    ops->MultiVecSetRandomValue((void **)multi_vec_y, num_vec_y, ops);
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
    ops->MultiVecPrint((void **)multi_vec_x, num_vec_x);
    printf ( "y\n" );
    ops->MultiVecPrint((void **)multi_vec_y, num_vec_y);


    double mvTmv[4];
    int    mv_s[2] = {0, 1};
    int    mv_e[2] = {2, 3};
    printf ( "MatDotMultiVec y = Ax\n" );
    printf ( "A\n" );
    CSR_PrintMat(A);
    ops->MatDotMultiVec((void *)A, (void **)multi_vec_x, (void **)multi_vec_y, 
	  mv_s, mv_e, ops);
    printf ( "y\n" );
    ops->MultiVecPrint((void **)multi_vec_y, num_vec_y);

    double alpha = 2.0, beta = -1.0;
    printf ( "MultiVecAxpby y = %0.4fx+%0.4fy\n", alpha, beta );
    ops->MultiVecAxpby(alpha, (void **)multi_vec_x, beta, (void **)multi_vec_y, 
	  mv_s, mv_e, ops);
    ops->MultiVecPrint((void **)multi_vec_y, num_vec_y);

    printf ( "MultiVecLinearComb y = x coord\n" );
    int ldcoord = 2;
    double coord[6] = {1.0, 1.0, -1.0, -1.0, 1.0, -1.0};
    printf ( "x\n" );
    ops->MultiVecPrint((void **)multi_vec_x, num_vec_x);
    printf ( "coord\n" );
    printf ( "%f\t%f\t%f\n", coord[0], coord[2], coord[4] );
    printf ( "%f\t%f\t%f\n", coord[1], coord[3], coord[5] );

    mv_s[0] = 0; mv_s[1] = 0;
    mv_e[0] = 2; mv_e[1] = 3;
    ops->MultiVecLinearComb((void **)multi_vec_x, (void **)multi_vec_y, 
	  mv_s, mv_e, coord, ldcoord, NULL, -1, ops);
    printf ( "y\n" );
    ops->MultiVecPrint((void **)multi_vec_y, num_vec_y);

    printf ( "MultiVecInnerProd value_ip\n" );
    double value_ip[6] = {1.0, 1.0, -1.0, -1.0, 1.0, -1.0};
    ops->MultiVecInnerProd((void **)multi_vec_x, (void **)multi_vec_y, 
	  (double *)&value_ip, "n", mv_s, mv_e, num_vec_x, ops);
    printf ( "value_ip\n" );
    printf ( "%f\t%f\t%f\n", value_ip[0], value_ip[2], value_ip[4] );
    printf ( "%f\t%f\t%f\n", value_ip[1], value_ip[3], value_ip[5] );

    printf ( "MultiVecSwap x<->y\n" );
    mv_s[0] = 0; mv_s[1] = 0;
    mv_e[0] = 2; mv_e[1] = 2;
    ops->MultiVecSwap((void **)multi_vec_x, (void **)multi_vec_y, 
	  mv_s, mv_e, ops);
    printf ( "x\n" );
    ops->MultiVecPrint((void **)multi_vec_x, num_vec_x);
    printf ( "y\n" );
    ops->MultiVecPrint((void **)multi_vec_y, num_vec_y);
    printf ( "Get and Restor x\n" );
    CSR_VEC *tmp;
    ops->GetVecFromMultiVec((void **)multi_vec_x, 1, (void *)&tmp);
    tmp->Entries[0] = 100.0;
    ops->RestoreVecForMultiVec((void **)multi_vec_x, 1, (void *)&tmp);
    printf ( "x\n" );
    ops->MultiVecPrint((void **)multi_vec_x, num_vec_x);



    //释放矩阵空间
    CSR_MatFree(&A);
    CSR_MatFree(&B);

    ops->MultiVecDestroy ((void ***)&multi_vec_x, num_vec_x, ops);
    ops->MultiVecDestroy ((void ***)&multi_vec_y, num_vec_y, ops);

    GCGE_OPS_Free(&ops);
//    mwTerm();
    return 0;
}
