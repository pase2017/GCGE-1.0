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

#if 1

void GCGE_Default_SolveLinearEquations(void *Matrix, void *b, void *x, GCGE_OPS *ops)
{
    //临时向量
    SLE_WS      *workspace = (SLE_WS*)(ops->solve_linear_equations_workspace);
    GCGE_INT    max_it     = workspace->cg_max_it;
    GCGE_DOUBLE rate       = workspace->cg_rate;
    void        **V_tmp    = workspace->vec_tmp;

    //临时向量
    void        *r, *p, *w;

    //CG迭代中用到的临时向量
    /* TODO 对于SLEPc中的多向量操作并不适用 */
    r = V_tmp[0]; //记录残差向量
    p = V_tmp[1];  //记录下降方向
    w = V_tmp[2];  //记录tmp=A*p

    GCGE_INT    niter = 0;
    GCGE_DOUBLE tmp1, tmp2, alpha, beta, rho1, rho2, error, last_error, pTw;

    ops->MatDotVec(Matrix,x,r); //tmp = A*x
    ops->VecAxpby(1.0, b, -1.0, r);  
    ops->VecInnerProd(r, r, &rho2);//用残量的模来判断误差
    error = sqrt(rho2);
    /* TODO */
    //这里应该判断以下如果error充分小就直接返回!!!
    printf("the initial residual: %e\n",error);
    //ops->VecAxpby(1.0, r, 0.0, p);
    niter = 1;
    last_error = error;
    printf("Rate = %e,  max_it = %d\n",rate,max_it);
    while((last_error/error >= rate)&&(niter<max_it))
    {
	if(niter == 1)
	{
	  //set the direction as r
	  ops->VecAxpby(1.0, r, 0.0, p);
	  //p = r;
	}
	else
	{
	  //printf("come here!\n");
	  //compute the value of beta
	  beta = rho2/rho1; 
	  //compute the new direction: p = r + beta * p
	  ops->VecAxpby(1.0, r, beta, p);
	}//end for if(niter == 1)
	//compute the vector w = A*p
	ops->MatDotVec(Matrix,p,w);
	//compute the value pTw = p^T * w 
	ops->VecInnerProd(p, w, &pTw);
	//compute the value of alpha
	//printf("pTw=%e\n",pTw);
	alpha = rho2/pTw; 
	//compute the new solution x = alpha * p + x
	ops->VecAxpby(alpha, p, 1.0, x);
	//compute the new residual: r = - alpha*w + r
	ops->VecAxpby(-alpha, w, 1.0, r);
	//set rho1 as rho2
	rho1 = rho2;
	//compute the new rho2
	ops->VecInnerProd(r, r, &rho2);  
	last_error = sqrt(rho2);      
//	printf("  niter= %d,  The current residual: %e\n",niter, last_error);
//	printf("  error=%e, last_error=%e,  last_error/error= %e\n",error, last_error, last_error/error);
	//update the iteration time
	niter++;   
    }
    printf("  niter= %d,  The residual: %e\n",niter, last_error);

}

#else
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

    printf("error: %e\n", error);
    printf ( "niter = %d\n", niter );
}
#endif

int main(int argc, char* argv[])
{
//    mwInit();
    srand((unsigned)time(NULL));

    //创建矩阵
    const char *file_A = "../data/A_5.txt";
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
    SLE_WS workspace = {1000, 1e-36, (void **)vec_tmp};
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
