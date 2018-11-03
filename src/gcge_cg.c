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
void GCGE_CG(void *Matrix, void *b, void *x, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
    //临时向量
    void        *r, *p, *w;
    GCGE_INT    niter = 0;
    GCGE_INT    max_it = para->cg_max_it, type = para->cg_type;
    GCGE_DOUBLE rate = para->cg_rate;
    GCGE_DOUBLE tmp1, tmp2, alpha, beta, rho1, rho2, error, last_error, pTw, wTw;
    void        **V_tmp = workspace->V_tmp;
    void        **RitzVec = workspace->RitzVec;

    //CG迭代中用到的临时向量
    ops->GetVecFromMultiVec(V_tmp, 1, &r); //记录残差向量
    ops->GetVecFromMultiVec(RitzVec, 0, &p); //记录下降方向
    ops->GetVecFromMultiVec(RitzVec, 1, &w); //记录tmp=A*p
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
        if(type == 1)
        {
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
        }
        else
        {    
            //这里我们把两个向量内积放在一起计算，这样可以减少一次的向量内积的全局通讯，提高可扩展性
            //compute the value pTw = p^T * w, wTw = w^T*w 
            ops->VecInnerProd(p, w, &pTw);
            ops->VecInnerProd(w, w, &wTw);

            alpha = rho2/pTw; 
            //compute the new solution x = alpha * p + x
            ops->VecAxpby(alpha, p, 1.0, x);
            //compute the new residual: r = - alpha*w + r
            ops->VecAxpby(-alpha, w, 1.0, r);
            //set rho1 as rho2
            rho1 = rho2;
            //compute the new rho2	
            rho2 = rho1 - 2.0*alpha*pTw + alpha*alpha*wTw;
        } 
        last_error = sqrt(rho2);      
        //printf("  niter= %d,  The current residual: %10.5f\n",niter, last_error);
        //printf("  error=%10.5f, last_error=%10.5f,  last_error/error= %10.5f\n",error, last_error, last_error/error);
        //update the iteration time
        niter++;   
    }//end while((last_error/error >= rate)&&(niter<max_it))
    ops->RestoreVecForMultiVec(V_tmp, 1, &r); 
    ops->RestoreVecForMultiVec(RitzVec, 0, &p);
    ops->RestoreVecForMultiVec(RitzVec, 1, &w);
}//end for the CG program




//统一进行块形式的CG迭代
//CG迭代求解 A * x =  RHS 
//W存储在 V(:,w_start:w_start+w_length), RHS存储在V_tmp(:,0:w_length)
// V_tmp = [r, p, w]的方式来组织
/**
 *  Matrix:    input,        A*x = RHS中的矩阵A
 *  RHS:       input,        前x_length个向量为要求解的x_length个方程组的右端项
 *  V:         input|output, 第x_start到x_start+x_length个向量为要求解的向量
 *                           input  为求解的初值
 *                           output 为得到的解
 *  x_start:   input,        V中解向量的起始位置
 *  x_length:  input,        要求解的方程组个数
 *  ops:       input,        操作
 *  para:      input,        参数
 *  workspace, input,        工作空间
 */
void GCGE_BCG(void *Matrix, void **RHS, void**V, GCGE_INT x_start, GCGE_INT x_length, 
          GCGE_OPS *ops, GCGE_PARA *para,GCGE_WORKSPACE *workspace)
{
  //临时向量
  void        *x,  *r, *p, *w, *b;
  void **CG_R = RHS, **CG_P = workspace->CG_p, **CG_W = workspace->evec, 
       **CG_X = V;
  GCGE_INT   id, idx, niter = 0;
  GCGE_INT    max_it = para->cg_max_it, type = para->cg_type;
  GCGE_DOUBLE rate = para->cg_rate;
  GCGE_DOUBLE tmp1, tmp2, alpha, beta;
  GCGE_DOUBLE *rho1 = workspace->subspace_dtmp, 
              *rho2 = rho1 + x_length, 
              *error = rho2 + x_length, 
              *last_error = error+x_length,
              *ptw = last_error + x_length,
              *rtr = ptw + x_length;
   GCGE_INT *unlock = workspace->subspace_itmp, num_unlock = 0; 


    for(idx=0; idx<x_length; idx++)
    {
      ops->GetVecFromMultiVec(CG_X, x_start+idx, &x);   //取出初始值 x      
      ops->GetVecFromMultiVec(CG_W, idx, &r); //取出暂存残差r的变量
      ops->MatDotVec(Matrix,x,r); //r = A*x
      ops->RestoreVecForMultiVec(CG_X, x_start+idx, &x);
      // b = b- r;
      ops->GetVecFromMultiVec(CG_R, idx, &b);   //取出右端项 b
      ops->VecAxpby(-1.0, r, 1.0, b);  
      ops->VecLocalInnerProd(b, b, rho2+idx);    //用残量的模来判断误差    
      ops->RestoreVecForMultiVec(CG_R, idx, &b); //把b相应的位置设为残差
      ops->RestoreVecForMultiVec(CG_W, idx, &r); //把r返回
      
      unlock[idx] = idx;
      num_unlock ++;
     }//算完初始的残差向量和残差
     //统一进行数据传输  
#if GCGE_USE_MPI
     MPI_Allreduce(MPI_IN_PLACE, rho2, x_length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
   for(idx=0;idx<x_length;idx++)
     error[idx] = sqrt(rho2[idx]);
    /**
     * ops->MatDotVec(Matrix,x,r); //tmp = A*x
    ops->VecAxpby(1.0, b, -1.0, r);  
    ops->VecInnerProd(r, r, &rho2);//用残量的模来判断误差
    error = sqrt(rho2); */
    /* TODO */
    //这里应该判断以下如果error充分小就直接返回!!!
    //printf("the initial residual: %e\n",error);
    //ops->VecAxpby(1.0, r, 0.0, p);
    niter = 1;
    //last_error = error; 
    //max_it = 40;
    //rate = 1e-5;
    //printf("Rate = %e,  max_it = %d\n",rate,max_it);
    while((num_unlock > 0)&&(niter < max_it))
    {
	if(niter == 1)
	{
	  for(id =0; id< num_unlock;id++)
	  {
	    idx = unlock[id];
	    ops->GetVecFromMultiVec(CG_P, idx, &p);   //取出初始值 x
	    ops->GetVecFromMultiVec(CG_R, idx, &r);   //取出右端项 r
	    //set the direction as r: p = r
	    ops->VecAxpby(1.0, r, 0.0, p);
	    ops->RestoreVecForMultiVec(CG_P, idx, &p);
	    ops->RestoreVecForMultiVec(CG_R, idx, &r); 
	  }//end for idx    
	}
	else
	{ 
	  for(id=0;id< num_unlock;id++)
	  {
	    idx = unlock[id];
	    //printf("come here!\n");
	    //compute the value of beta
	    beta = rho2[idx]/rho1[idx]; 
	    //compute the new direction: p = r + beta * p
	    ops->GetVecFromMultiVec(CG_R, idx, &r);   //取出右端项 r
	    ops->GetVecFromMultiVec(CG_P, idx, &p);   //取出初始值 x
	    //p = r + beta * p
	    ops->VecAxpby(1.0, r, beta, p);
	    ops->RestoreVecForMultiVec(CG_P, idx, &p);
	    ops->RestoreVecForMultiVec(CG_R, idx, &r);
	  }//end for idx
	}//end for if(niter == 1)
	//compute the vector w = A*p和p^Tw = ptw
	for(id = 0; id < num_unlock; id++)
	{
	  idx = unlock[id];
	  ops->GetVecFromMultiVec(CG_P, idx, &p);   //取出 p
	  ops->GetVecFromMultiVec(CG_W, idx, &w);   //取出 w
	  ops->MatDotVec(Matrix,p,w);  //w = A*p
	  //做局部内积（在每个进程内部做局部的内积）
	  ops->VecLocalInnerProd(p,w,ptw+idx);
	  //ops->VecLocalInnerProd(w,w,wtw+idx);	   
	  ops->RestoreVecForMultiVec(CG_P, idx, &p); 
	  ops->RestoreVecForMultiVec(CG_W, idx, &w);
	}//end for idx
	//这里我们把两个向量内积放在一起计算，这样可以减少一次的向量内积的全局通讯，提高可扩展性
	//compute the value pTw = p^T * w
	//应该要写一个专门的子程序来做这个统一的内积
	//与张宁商量一下如何实现
#if GCGE_USE_MPI
     MPI_Allreduce(MPI_IN_PLACE, ptw, x_length,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif     
     //现在来计算: alpha = rho2/ptw, x = alpha*p +x , r = -alpha*w + r, 并且同时计算r的局部内积
     for(id=0;id<num_unlock;id++)
     { 
       idx = unlock[id];
       //计算alpha
       alpha = rho2[idx]/ptw[idx];
       //compute the new solution x = alpha * p + x
       ops->GetVecFromMultiVec(CG_P, idx, &p);   //取出 p

       ops->GetVecFromMultiVec(CG_X, x_start+idx, &x); //取出初始值 x 
       //  x = alpha*p +x 
       ops->VecAxpby(alpha, p, 1.0, x);
       ops->RestoreVecForMultiVec(CG_P, idx, &p); 
       ops->RestoreVecForMultiVec(CG_X, x_start+idx, &x);
       
       
       ops->GetVecFromMultiVec(CG_R, idx, &r);   //取出 r
       ops->GetVecFromMultiVec(CG_W, idx, &w);   //取出 w
       //r = -alpha*w + r
       ops->VecAxpby(-alpha, w, 1.0, r);
       //set rho1 as rho2
       rho1[idx] = rho2[idx];
       //计算新的r的局部内积
       ops->VecLocalInnerProd(r, r, rho2+idx);       
       ops->RestoreVecForMultiVec(CG_R, idx, &r); 
       ops->RestoreVecForMultiVec(CG_W, idx, &w);

     }//end for idx
     //统一进行数据传输  
#if GCGE_USE_MPI
     MPI_Allreduce(MPI_IN_PLACE, rho2, x_length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
   for(idx=0;idx<x_length;idx++)
     last_error[idx] = sqrt(fabs(rho2[idx]));  
   
   //printf("niter= %d, Error[%d]=%e, rho2[%d]=%e \n", niter, idx, last_error[idx], idx, rho2[idx]);
    //下面进行收敛性的判断
     num_unlock = 0;
     for(idx=0;idx<x_length;idx++)
     { 
       if(last_error[idx]/error[idx] >= rate)
       {
	 unlock[idx] = idx;
	 num_unlock ++;
	}//end if(last_error[idx]/error[idx] >= rate)       
     }//end for idx
    //printf("  niter= %d,  The current residual: %10.5f\n",niter, last_error);
    //printf("  error=%10.5f, last_error=%10.5f,  last_error/error= %10.5f\n",error, 
    //last_error, last_error/error);
    //update the iteration time
    niter++;       
    //printf("\n");   
    }//end while((last_error/error >= rate)&&(niter<max_it))
}//end for this subprogram
