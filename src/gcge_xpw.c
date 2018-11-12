/*
 * =====================================================================================
 *
 *       Filename:  gcge_xpw.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2018年09月24日 09时50分16秒
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

#include "gcge_xpw.h"


/**
 * @brief  计算 X = V * a 
 *
 *   n_V = workspace->dim_xpw;
 *   n_X = workspace->dim_xpw;
 *   lda = n_V;
 *
 *   start = { 0, n_V }
 *   end   = { 0, n_X }
 *
 *   ops->MultiVecLinearComb(V, X, start, end, a, lda, NULL, -1, ops);
 *
 * @param V
 * @param subspace_evec
 * @param RitzVec
 * @param ops
 * @param workspace
 */
//计算RitzVec = V*subspace_evec
//V是用来计算Ritz向量的基底向量组(GCGE_Vec),RitzVec是计算得到的Ritz向量(GCGE_Vec)
void GCGE_ComputeX(void **V, GCGE_DOUBLE *subspace_evec, void **X, 
        GCGE_INT size_X, GCGE_INT size_V,
        GCGE_OPS *ops, GCGE_WORKSPACE *workspace)
{
    /* TODO */
#if 0
    GCGE_DOUBLE    t1 = GCGE_GetTime();
#endif
    GCGE_INT i = 0;
    GCGE_INT dim_v = size_V;
    GCGE_INT dim_x = size_X;
    GCGE_INT start[2];
    GCGE_INT end[2];

    //使用V(:,start_x:dim_v)线性组合得到X(:,start_x:dim_x)
    //线性组合系数从start行start列开始
    GCGE_INT start_x = workspace->num_soft_locked_in_last_iter;
    start[0] = start_x;
    end[0]   = dim_v;
    start[1] = start_x;
    end[1]   = dim_x;
    ops->MultiVecLinearComb(V, X, start, end, subspace_evec+start_x*dim_v+start_x,
            dim_v, NULL, 0, ops);
#if 0 
    GCGE_DOUBLE t2 = GCGE_GetTime();
    //统计计算Rtiz向量时间
    solver->para->stat_para->PartTimeTotal->GetX_Time += t2-t1;
    solver->para->stat_para->PartTimeInOneIte->GetX_Time = t2-t1;
#endif
}

/**
 * @brief 计算W：求解 A * W = lambda * B * X
 *        只对 X 中编号在 unlock 中的向量进行此操作
 *  
 *  形式参数：A, B, eval, V, w_start, w_length, unlock, void *V_tmp
 *  CG 迭代需要3个向量的临时空间,存储rhs使用
 *  临时空间V_tmp中需要最少4个向量的空间
 *  TODO
 *  1. 之后CG会改成与ops结合
 *  2. CG需要的参数是否需要workspace
 *
 *  for idx = 0 : w_length
 *      x_current = unlock[idx]
 *      w_current = w_start + idx
 *
 *      计算 rhs = B * X[x_current]
 *      if B == NULL
 *          VecAxpby(1.0, X[x_current], 0.0, rhs);
 *      else
 *          MatDotVec(B, X[x_current], rhs);
 *      end
 *
 *      计算 rhs = lambda * B * X[x_current]
 *      VecAxpby(0.0, rhs, lambda, rhs);
 *
 *      计算 W[w_current] : A * W[w_current] = rhs
        GCGE_CG(A, rhs, W[w_current], ops, para, V_tmp);
 *
 *  end
 *
 * @param A
 * @param B
 * @param V
 * @param eval
 * @param ops
 * @param para
 * @param workspace
 */
//用未收敛的特征向量作为X求解线性方程组 A W = \lambda B X
//A,B表示矩阵A,B(GCGE_MATRIX),V表示向量组(GCGE_Vec)
void GCGE_ComputeW(void *A, void *B, void **V, GCGE_DOUBLE *eval, 
        GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
    GCGE_INT w_start  = workspace->dim_xp;     //W向量在V空间中的起始位置
    GCGE_INT w_length = workspace->unconv_bs;  //W的向量个数
    GCGE_INT *unlock  = workspace->unlock;     //未收敛的编号
    void     **V_tmp  = workspace->V_tmp;      //临时空间，用于存储rhs以及CG迭代使用
    void     *rhs;                             //线性方程组求解时的右端项
    void     *x_current;                       //当前操作的X向量
    void     *w_current;                       //当前操作的W向量


    GCGE_INT idx = 0;
    GCGE_INT x_current_idx = 0;
    GCGE_INT w_current_idx = 0;
    //先把W方程的右端项，初值都计算好
    for(idx=0; idx<w_length; idx++)
    {
        //rhs = V_tmp(:,0);
        ops->GetVecFromMultiVec(V_tmp, idx, &rhs);
        x_current_idx = unlock[idx];
        w_current_idx = w_start + idx;	
        ops->GetVecFromMultiVec(V, x_current_idx, &x_current);
        ops->GetVecFromMultiVec(V, w_current_idx, &w_current);
        // B == NULL，那么 rhs = lambda * x_current
        // B != NULL，那么 rhs = lambda * B * x_current
        if(B == NULL)
        {
            ops->VecAxpby(1.0, x_current, 0.0, rhs);
        }
        else
        {
            ops->MatDotVec(B, x_current, rhs);
        }
        //对lambda的用途：
        //赋初值
        if(fabs(eval[x_current_idx]) <= 1.0)
        {
            //如果lambda <= 1.0, rhs = lambda * rhs
            ops->VecAxpby(0.0, rhs, eval[x_current_idx], rhs);
            //给CG迭代赋初值： w = x; (A*x = lambda *B *x)
            ops->VecAxpby(1.0, x_current, 0.0, w_current);	  
        }
        else
        {
            //当lambda>1.0的时候，定义初值解为: 1/lambda*x; 
            //线性方程: (A*(x/lambda) = B * x)
            ops->VecAxpby(1.0/eval[x_current_idx], x_current, 0.0, w_current);	  
        }
        //再把x_current 和 w_current 返回回去
        ops->RestoreVecForMultiVec(V, x_current_idx, &x_current);
        ops->RestoreVecForMultiVec(V, w_current_idx, &w_current);
        //把rhs返回回去
        ops->RestoreVecForMultiVec(V_tmp, idx, &rhs);
    }//for(idx=0; idx<w_length; idx++)

    //下面统一进行线性方程的求解           
    //用户提供了解法器，就用用户的
    //W存储在 V(:,w_start:w_start+w_length), RHS存储在V_tmp(:,0:w_length)
    if(ops->LinearSolver)
    {
        //这里是直接调用用户的解法器，期望未来用户也能提供块形式的解法器（统一进行数据广播的形式）
        for(idx=0; idx<w_length; idx++)
        {
            ops->GetVecFromMultiVec(V_tmp, idx, &rhs);
            w_current_idx = w_start + idx;	
            ops->GetVecFromMultiVec(V, w_current_idx, &w_current);
            ops->LinearSolver(A, rhs, w_current, ops);
            ops->RestoreVecForMultiVec(V, w_current_idx, &w_current);
            ops->RestoreVecForMultiVec(V_tmp, idx, &rhs);
        }//end for idx_p      
	    //如果有外部求解器, 解完后继续用BCG求解. 如果不要用BCG, 可以将cg_max_it设为0
        //GCGE_BCG(A, V_tmp, V, w_start,w_length, ops, para,workspace);
    }
    else
    {
        //统一进行块形式的CG迭代
        //CG迭代求解 A * W =  RHS 
        //W存储在 V(:,w_start:w_start+w_length), RHS存储在V_tmp(:,0:w_length)
        if(para->if_use_bcg == 1)
        {
            GCGE_BCG(A, V_tmp, V, w_start,w_length, ops, para,workspace);
        }
        else
        {
            //如果不使用BCG, 因为在GCGE_CG中还要取出三个向量使用, 所以使用slepc时不能使用这种普通CG
            for(idx=0; idx<w_length; idx++)
            {
                ops->GetVecFromMultiVec(V_tmp, idx, &rhs);
                w_current_idx = w_start + idx;	
                ops->GetVecFromMultiVec(V, w_current_idx, &w_current);
                GCGE_CG(A, rhs, w_current, ops, para, workspace);
                ops->RestoreVecForMultiVec(V, w_current_idx, &w_current);
                ops->RestoreVecForMultiVec(V_tmp, idx, &rhs);
            }//end for idx_p      
        }
    }//end for LinearSolver    
}//end for this subprogram


/**
 * @brief 计算P的子空间系数，并进行小规模正交化
 * 符号解释: XX: 表示新的X向量对现有的基向量组中X的依赖关系，
 *           PX： 表示新的P向量对现有的基向量组中X的依赖关系
 *			其他命名类同
 *        再线性组合得到P
 *        XX  O  拷贝  XX  O   正交化  XX  PX  线性组合         PX
 *        XP  O  ===>  XP  XP  =====>  XP  PP  =======> P = V * PP
 *        XW  O        XW  XW          XW  PW                   PW
 *
 *  形式参数：coef, ldc, unlock, p_start, p_ncols, 
 *            xx_nrows: 基向量中X向量的个数,
 *            xx_ncols: 更新后基向量中X的个数,由CheckMultiplity 函数计算得到
 *            V, V_tmp: 线性组合的临时空间，需要p_ncols的向量空间,至少block_size,
 *            ops
 *
 * >coef :
 *
 *  XX  O      XX  O
 *  XP  O  =>  XP  XP
 *  XW  O      XW  XW
 *
 * >赋0
 *  PX的指针 px = coef + ldc * p_start
 *  for idx=0; idx<p_ncols; idx++
 *      memset(px+idx*ldc, 0.0, xx_nrows*sizeof(double))
 *  end
 *
 * >复制
 *  PP的指针 pp = coef + ldc * p_start + xx_nrows
 *  XP的指针 xp = coef + xx_nrows
 *  for idx_p=0; idx_p<p_ncols; idx_p++
 *      current_x = unlock[idx_p]
 *      memcpy(pp+idx_p*ldc, xp+current_x*ldc, (ldc-xx_nrows)*sizeof(double))
 *  end
 *
 * >对coef的p_start到p_end进行正交化
 *  p_end = p_start + p_ncols
 *  GCGE_OrthogonalSubspace(coef, ldc, p_start, &p_end, NULL, -1, orth_para)
 *  p_ncols = p_end - p_start
 *
 * >V线性组合得到X中看作重特征值的特征向量，存放在V_tmp中
 *  x_end_in_V_tmp = p_start-nev
 *  mv_s = { 0, 0 }
 *  mv_e = { ldc, x_end_in_V_tmp }
 *  MultiVecLinearComb(V, V_tmp, mv_s, mv_e, px, ldc, NULL, -1, ops)
 * >V线性组合得到P向量，存放在V_tmp中
 *  p_end_in_V_tmp  = x_end_in_V_tmp + p_ncols
 *  mv_s = { 0, x_end_in_V_tmp }
 *  mv_e = { ldc, p_end_in_V_tmp }
 *  MultiVecLinearComb(V, V_tmp, mv_s, mv_e, px, ldc, NULL, -1, ops)
 *  把V_tmp与V中X中看作重特征值的特征向量的部分进行指针交换
 *  把V_tmp与V中P的部分进行指针交换
 *  mv_s = { nev, 0 }
 *  mv_e = { p_end, p_end_in_V_tmp }
 *  MultiVecSwap(V, V_tmp, mv_s, mv_e, ops)
 *  mv_s = { nev, 0 }
 *  mv_e = { p_start, p_ncols }
 *  MultiVecSwap(V, V_tmp, mv_s, mv_e, ops)
 *
 * > out : 更新 dim_xp
 *         更新 V
 *
 */
//获取P向量组
//V表示用于计算的基底向量组(GCGE_Vec)
void GCGE_ComputeP(GCGE_DOUBLE *subspace_evec, void **V, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
    //t1,t2用于统计时间
    GCGE_DOUBLE t1 = 0.0;
    GCGE_DOUBLE t2 = 0.0;
#if GET_PART_TIME
    t1 = GCGE_GetTime();
#endif
    //要给的形式参数
    GCGE_DOUBLE *coef    = subspace_evec; //线性组合系数
    GCGE_INT    ldc      = workspace->dim_xpw; //系数矩阵ldc的leading dimension
    GCGE_INT    *unlock  = workspace->unlock; //未收敛的编号,要对X中unlock的列进行操作
    GCGE_INT    p_start  = workspace->dim_x; //P向量(在V中)的起始位置
    GCGE_INT    p_ncols  = workspace->unconv_bs; //P向量的个数
    GCGE_INT    xx_ncols = workspace->dim_x; //更新后基向量中X的个数,由CheckMultiplity 函数计算得到
    GCGE_INT    xx_nrows = workspace->last_dim_x; //基向量中, 现有的X的向量个数,
    void        **V_tmp  = workspace->V_tmp;
    GCGE_INT    mv_s[2]  = { 0,0 };
    GCGE_INT    mv_e[2]  = { 0,0 };

    //函数中用到的参数定义
    GCGE_INT    current_p = 0; //当前P向量的位置
    GCGE_INT    current_x = 0; //当前X向量的位置
    GCGE_INT    p_end     = 0; //P向量的终点位置

    /* XX  PX 
     * XP  PP 
     * XW  PW 
     * px表示上述矩阵中PX的起始位置
     * pp表示上述矩阵中PP的起始位置
     * xp表示上述矩阵中XP的起始位置 */
     //px: 表示P中第一个向量对X的依赖系数在subspace_evec中的起始位置
    GCGE_DOUBLE *px = coef + p_start * ldc;
    //pp: 表示P中第一个向量对P的依赖系数在subspace_evec中的起始位置
    GCGE_DOUBLE *pp = px   + xx_nrows;
    //xp: 表示X中第一个向量对P的依赖系数在subspace_evec中的起始位置
    GCGE_DOUBLE *xp = coef + xx_nrows;

    //初始化为0并进行复制,详细见函数名前的注释
    memset(px, 0.0, p_ncols*ldc*sizeof(GCGE_DOUBLE));
    for(current_p=0; current_p<p_ncols; current_p++)
    {
        current_x = unlock[current_p];
        memcpy(pp+current_p*ldc, xp+current_x*ldc, (ldc-xx_nrows)*sizeof(GCGE_DOUBLE));
    }
    //现在得到了与[X,P]相对应的coeff矩阵的初始向量
    //对coef的p_start到p_end列进行正交化，会修改p_end
    /* 对稠密矩阵 C 表示连续的已收敛的locked部分
     * 此时子空间特征向量矩阵有下面的形式
     *
     *  CC  XC  PC
     *  CN  XN  PN
     *  CP  XP  PP
     *  CW  XW  PW
     *
     * XC 表示 XX 中 未收敛的特征向量的连续的已收敛的locked部分
     * PC 表示 PX 中 连续的已收敛的locked部分
     * XN 表示 [XX^T,XP^T,XW^T]^T 中 
     *         未收敛的特征向量的第一个未收敛的及之后的部分
     * PN 表示 PX 中 第一个未收敛的及之后的部分
     * CC,CN,CP,CW表示 [XX^T,XP^T,XW^T]^T 中已收敛的特征向量的部分
     *
     * 对其中的
     *
     *  XN  PN
     *  XP  PP
     *  XW  PW
     *
     * 部分进行正交化
     *
     */
    // orth_start_row 表示从稠密矩阵的哪一行开始正交化,即PN的起始行号
    GCGE_INT orth_start = workspace->num_soft_locked_in_last_iter;
    // orth_length 表示 进行正交化的行数
    GCGE_INT orth_length = ldc - orth_start;
    p_start -= orth_start;
    p_end = p_start + p_ncols;
    if(strcmp(para->p_orth_type, "scbgs") == 0)
    {
        GCGE_SCBOrthogonalSubspace(coef + orth_start*ldc + orth_start, 
                ldc, orth_length, p_start, &p_end, NULL, -1, 
                para->orth_para, workspace, ops);
    }
    else if(strcmp(para->p_orth_type, "bgs") == 0)
    {
        GCGE_BOrthogonalSubspace(coef + orth_start*ldc + orth_start, 
                ldc, orth_length, p_start, &p_end, NULL, -1, 
                para->orth_para, workspace->subspace_dtmp, ops);
    }
    else
    {
        GCGE_OrthogonalSubspace(coef + orth_start*ldc + orth_start, ldc, 
                orth_length, p_start, &p_end, NULL, -1, para->orth_para);
    }
    p_ncols = p_end - p_start;
    p_end += orth_start;

#if GET_PART_TIME
    t2 = GCGE_GetTime();
    para->stat_para->part_time_one_iter->p_orth_time = t2-t1;
    para->stat_para->part_time_total->p_orth_time += t2-t1;
#endif
    //V线性组合得到P向量，存放在V_tmp中
    //V线性组合得到X中看作重特征值的部分特征向量，存放在V_tmp中
    GCGE_INT nev = para->nev;
    GCGE_INT p_end_in_V_tmp = 0;
    p_end_in_V_tmp = p_end - nev;

    //线性组合的基底向量为V(:,orth_start:ldc)
    mv_s[0] = orth_start;
    mv_e[0] = ldc;
    mv_s[1] = 0;
    mv_e[1] = p_end_in_V_tmp;
    //线性组合的系数从第orth_start行开始
    ops->MultiVecLinearComb(V, V_tmp, mv_s, mv_e, coef+nev*ldc+orth_start, 
            ldc, NULL, -1, ops);

    //把V_tmp与V中P的部分进行指针交换
    //把V_tmp与V中X看作重特征值的部分特征向量的部分进行指针交换
    mv_s[0] = nev;
    mv_s[1] = 0;
    mv_e[0] = p_end;
    mv_e[1] = p_end_in_V_tmp;
    ops->MultiVecSwap(V, V_tmp, mv_s, mv_e, ops);

    //更新workspace中的dim_xp
    workspace->dim_xp = p_end;
#if GET_PART_TIME
    t1 = GCGE_GetTime();
    para->stat_para->part_time_one_iter->p_axpy_time = t1-t2;
    para->stat_para->part_time_total->p_axpy_time += t1-t2;
#endif

}

