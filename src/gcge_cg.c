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
    void        *r, *p, *tmp;
    GCGE_INT    niter = 0;
    GCGE_INT    max_it = para->cg_max_it;
    GCGE_DOUBLE rate = para->cg_rate;
    GCGE_DOUBLE tmp1, tmp2, alpha, beta, error, last_error;

    //CG迭代中用到的临时向量
    r   = V_tmp[1];
    p   = V_tmp[2];
    tmp = V_tmp[3];
  
    ops->MatDotVec(Matrix,x,p); //tmp1 = A*x0
 
    ops->VecAxpby(1.0, b, 0.0, r);
    ops->VecAxpby(-1.0, p, 1.0, r);//r=b-p=r-p
    ops->VecInnerProd(r, r, &error);//用残量的模来判断误差
    error = sqrt(error);
    ops->VecAxpby(1.0, r, 0.0, p);

    do{
        ops->MatDotVec(Matrix,p,tmp);//tmp = A*p
        ops->VecInnerProd(r,p,&tmp1);
        ops->VecInnerProd(p,tmp,&tmp2);    
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
    
        //beta = -(r,tmp)/(p,tmp)
        ops->VecInnerProd(r,tmp,&tmp1);
        ops->VecInnerProd(p,tmp,&tmp2);    
        beta = -tmp1/tmp2;
        ops->VecAxpby(1.0, r, beta, p);    

        niter++; 
    }while((error/last_error >= rate)&&(niter<max_it));

}

