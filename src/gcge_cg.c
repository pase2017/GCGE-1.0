/*
 * =====================================================================================
 *
 *       Filename:  gcge_cg.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2018年09月24日 09时57分13秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#include <math.h>

#include "gcge_cg.h"

//CG迭代求解Matrix * x = b
//Matrix是用于求解的矩阵(GCGE_MATRIX),b是右端项向量(GCGE_VEC),x是解向量(GCGE_VEC)
//用时，V_tmp+1
void GCGE_CG(void *Matrix, void *b, void *x, GCGE_OPS *ops, GCGE_PARA *para, void **V_tmp)
{
    //临时向量
    void        *r, *p, *w;
    GCGE_INT    niter = 0;
    GCGE_INT    max_it = para->cg_max_it;
    GCGE_DOUBLE rate = para->cg_rate;
    GCGE_DOUBLE tmp1, tmp2, alpha, beta, rho1, rho2, error, last_error, pTw;

    //CG迭代中用到的临时向量
    r = V_tmp[1]; //记录残差向量
    p = V_tmp[2];  //记录下降方向
    w = V_tmp[3];  //记录tmp=A*p
    ops->MatDotVec(Matrix,x,r); //tmp = A*x
    ops->VecAxpby(1.0, b, -1.0, r);  
    ops->VecInnerProd(r, r, &rho2);//用残量的模来判断误差
    error = sqrt(rho2);
    /* TODO */
    //这里应该判断以下如果error充分小就直接返回!!!
    //printf("the initial residual: %e\n",error);
    //ops->VecAxpby(1.0, r, 0.0, p);
    niter = 1;
    last_error = error;
    //max_it = 40;
    //rate = 1e-5;
    //printf("Rate = %e,  max_it = %d\n",rate,max_it);
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
        //printf("  niter= %d,  The current residual: %10.5f\n",niter, last_error);
        //printf("  error=%10.5f, last_error=%10.5f,  last_error/error= %10.5f\n",error, last_error, last_error/error);
        //update the iteration time
        niter++;   
    }

}

