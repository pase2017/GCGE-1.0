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

typedef struct SLE_WS_ {

   GCGE_INT    cg_max_it;
   GCGE_DOUBLE cg_rate;
   void**      vec_tmp; /* 长度为3的多向量指针 */

}SLE_WS;

void GCGE_Default_SolveLinearEquations(void *mat, void *b, void *x, struct GCGE_OPS_ *ops)
{
    //临时向量
    SLE_WS      *workspace = (SLE_WS*)(ops->solve_linear_equations_workspace);
    GCGE_INT    max_it     = workspace->cg_max_it;
    GCGE_DOUBLE rate       = workspace->cg_rate;
    void        **V_tmp    = workspace->vec_tmp;

    void        *r, *p, *tmp;
    GCGE_INT    niter = 0;
    GCGE_DOUBLE tmp1, tmp2, alpha, beta, error, last_error;

    //CG迭代中用到的临时向量
    /* TODO 对于SLEPc中的多向量操作并不适用 */
    r   = V_tmp[0];
    p   = V_tmp[1];
    tmp = V_tmp[2];
  
    ops->MatDotVec(mat,x,p); //tmp1 = A*x0
 
    ops->VecAxpby(1.0, b, 0.0, r);
    ops->VecAxpby(-1.0, p, 1.0, r);//r=b-p=r-p
    ops->VecInnerProd(r, r, &error);//用残量的模来判断误差
    error = sqrt(error);
    ops->VecAxpby(1.0, r, 0.0, p);

    /* TODO 需要判断error是否已经很小, 否则可能出现计算alpha为nan */
    do{
        ops->MatDotVec(mat,p,tmp);//tmp = A*p
        ops->VecInnerProd(r,p,&tmp1);
        ops->VecInnerProd(p,tmp,&tmp2);    
	/* TODO 需要判断分母是否接近于0 */
        alpha = tmp1/tmp2;
        ops->VecAxpby(alpha, p, 1.0, x);   
        ops->VecAxpby(-alpha, tmp, 1.0, r);
        last_error = error;
        ops->VecInnerProd(r, r, &error);   //用残量的模来判断误差
        error = sqrt(error);
    
        if(error/last_error < rate)
        {
            //printf("error: %e, last_error: %e, error/last_error: %e, rate: %e\n", 
            //    error, last_error, error/last_error, rate);
            break;
        }
    
        /*beta = -(r,tmp)/(p,tmp)*/
        ops->VecInnerProd(r,tmp,&tmp1);
        ops->VecInnerProd(p,tmp,&tmp2);    
        beta = -tmp1/tmp2;
        ops->VecAxpby(1.0, r, beta, p);    

        niter++; 
    }while((error/last_error >= rate)&&(niter<max_it));

    printf ( "niter = %d\n", niter );
}

int main(int argc, char* argv[])
{
//    mwInit();
    srand((unsigned)time(NULL));

    //创建矩阵
    const char *file_A = "../data/A_3.txt";
//    const char *file_A = "../data/testA";
    CSR_MAT *A = CSR_ReadMatFile(file_A);
    CSR_VEC *rhs, *x, **vec_tmp;

    GCGE_OPS *ops;
    GCGE_OPS_Create(&ops);
    GCGE_CSR_SetOps(ops);
    GCGE_OPS_Setup(ops);

    ops->BuildVecByMat((void *)A, (void **)&rhs);
    ops->BuildVecByMat((void *)A, (void **)&x);
    ops->BuildMultiVecByMat((void *)A, (void ***)&vec_tmp, 3, ops);

    ops->SolveLinearEquations = GCGE_Default_SolveLinearEquations;
    SLE_WS workspace = {30, 1e-2, (void **)vec_tmp};
    ops->solve_linear_equations_workspace = &workspace;

    double norm;
    for (int i = 0; i < 10; ++i)
    {
       ops->VecSetRandomValue((void *)rhs);
       ops->VecSetRandomValue((void *)x);
       ops->SolveLinearEquations((void *)A, (void *)rhs, (void *)x, ops);

       ops->MatDotVec((void *)A, (void *)x, (void *)vec_tmp[0]);
       ops->VecAxpby(-1, (void *)vec_tmp[0], 1, (void *)rhs);
       ops->VecNorm(rhs, &norm, ops);

       printf ( "norm of res = %e\n", norm );
    }


    //释放矩阵空间
    CSR_MatFree(&A);
    ops->FreeMultiVec((void ***)&vec_tmp, 3, ops);
    ops->FreeVec((void **)&rhs);
    ops->FreeVec((void **)&x);

    GCGE_OPS_Free(&ops);

//    mwTerm();
    return 0;
}
