/*
 * =====================================================================================
 *
 *       Filename:  test_cg.c
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

    CSR_PrintMat(A);
    CSR_PrintMat(B);

    CSR_VEC *rhs;
    CSR_BuildVecByMat(A, &rhs);

    srand( (unsigned)time( NULL ) ); 
    int i;
    for(i=0; i<rhs->size; i++)
    {
       rhs->Entries[i] = 10*rand()/(RAND_MAX+1.0);
    }

    CSR_VEC *x;
    CSR_BuildVecByMat(A, &x);

    CSR_VEC **vec_tmp = (CSR_VEC**)malloc(4*sizeof(CSR_VEC*));
    for(i=0; i<4; i++)
    {
        CSR_BuildVecByMat(A, vec_tmp+i);
    }

    printf ( "Ax=rhs\n" );
    printf ( "A\n" );
    CSR_PrintMat(A);
    printf ( "rhs\n" );
    CSR_PrintVec(rhs);

    GCGE_OPS *ops;
    GCGE_OPS_Create(&ops);
    GCGE_CSR_SetOps(ops);
    GCGE_PARA *para;
    GCGE_PARA_Create(&para);
    GCGE_PARA_Setup(para);
    /* Should op->SolveLinearEquations */
    /* Should ls_para linear_solver_para */
    GCGE_CG((void *)A, (void *)rhs, (void *)x, ops, para, (void **)vec_tmp);

    printf ( "x\n" );
    CSR_PrintVec(x);
    printf ( "rhs - Ax\n" );
    ops->MatDotVec((void *)A, (void *)x, (void *)vec_tmp[0]);
    ops->VecAxpby(-1, (void *)vec_tmp[0], 1, (void *)rhs);
    CSR_PrintVec(rhs);

    GCGE_PARA_Free(&para);
    GCGE_OPS_Free(&ops);

    //释放矩阵空间
    CSR_MatFree(&A);
    CSR_MatFree(&B);
    CSR_VecFree(&rhs); rhs = NULL;
    CSR_VecFree(&x);   x   = NULL;

    for(i=0; i<4; i++)
    {
        CSR_VecFree(vec_tmp+i);
    }
    free(vec_tmp); vec_tmp = NULL;

//    mwTerm();
    return 0;
}
