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
void GCGE_ComputeX(void **V, GCGE_DOUBLE *subspace_evec, void **RitzVec, 
        GCGE_OPS *ops, GCGE_WORKSPACE *workspace)
{
    /* TODO */
#if 0
    GCGE_DOUBLE    t1 = GCGE_GetTime();
#endif
    GCGE_INT i = 0;
    GCGE_INT dim_v = workspace->dim_xpw;
    GCGE_INT dim_x = workspace->dim_x;
    GCGE_INT start[2];
    GCGE_INT end[2];
    start[0] = 0;
    end[0]   = dim_v;
    start[1] = 0;
    end[1]   = dim_x;
    ops->MultiVecLinearComb(V, RitzVec, start, end, subspace_evec, dim_v, NULL, 0, ops);
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
    GCGE_INT w_start  = workspace->dim_xp;
    GCGE_INT w_length = workspace->unconv_bs;
    GCGE_INT *unlock  = workspace->unlock;
    void     **V_tmp  = workspace->V_tmp;
    void     *rhs;
    void     *x_current;
    void     *w_current;
    ops->GetVecFromMultiVec(V_tmp, 0, &rhs);

    GCGE_INT idx = 0;
    GCGE_INT x_current_idx = 0;
    GCGE_INT w_current_idx = 0;

    for(idx=0; idx<w_length; idx++)
    {
        x_current_idx = unlock[idx];
        w_current_idx = w_start + idx;
        ops->GetVecFromMultiVec(V, x_current_idx, &x_current);
        ops->GetVecFromMultiVec(V, w_current_idx, &w_current);
        if(B == NULL)
        {
            ops->VecAxpby(1.0, x_current, 0.0, rhs);
        }
        else
        {
            ops->MatDotVec(B, x_current, rhs);
        }
        ops->VecAxpby(0.0, rhs, eval[x_current_idx], rhs);
        GCGE_CG(A, rhs, w_current, ops, para, V_tmp);
        ops->RestoreVecForMultiVec(V, x_current_idx, &x_current);
        ops->RestoreVecForMultiVec(V, w_current_idx, &w_current);
    }
    ops->RestoreVecForMultiVec(V_tmp, 0, &rhs);
}


/**
 * @brief 计算P的子空间系数，并进行小规模正交化
 *        再线性组合得到P
 *
 *  形式参数：coef, ldc, unlock, p_start, p_ncols, 
 *            xx_nrows: 基向量中X向量的个数,
 *            xx_ncols: 更新后基向量中X的个数,由CheckMultiplity 函数计算得到
 *            V, V_tmp: 线性组合的临时空间，需要p_ncols的向量空间,至少block_size,
 *            ops
 *
 * coef :
 *
 *  XX  PX      XX  O
 *  XP  PP  =>  XP  XP
 *  XW  PW      XW  XW
 *
 *  PX的指针 px = coef + ldc * p_start
 *  赋0
 *  for idx=0; idx<p_ncols; idx++
 *      memset(px+idx*ldc, 0.0, xx_nrows*sizeof(double))
 *  end
 *  复制
 *  PP的指针 pp = coef + ldc * p_start + xx_nrows
 *  XP的指针 xp = coef + xx_nrows
 *  for idx_p=0; idx_p<p_ncols; idx_p++
 *      current_x = unlock[idx_p]
 *      memcpy(pp+idx_p*ldc, xp+current_x*ldc, (ldc-xx_nrows)*sizeof(double))
 *  end
 *  p_end = p_start + p_ncols
 *  对coef的p_start到p_end进行正交化
 *  GCGE_OrthogonalSubspace(coef, ldc, p_start, &p_end, NULL, -1, orth_para)
 *  p_ncols = p_end - p_start
 *  V线性组合得到P向量，存放在V_tmp中
 *  mv_s = { 0, 0 }
 *  mv_e = { ldc, p_ncols }
 *  MultiVecLinearComb(V, V_tmp, mv_s, mv_e, px, ldc, NULL, -1, ops)
 *  把V_tmp与V中P的部分进行指针交换
 *  mv_s = { p_start, 0 }
 *  mv_e = { p_end, p_ncols }
 *  MultiVecSwap(V, V_tmp, mv_s, mv_e, ops)
 *
 * @param subspace_evec
 * @param V
 * @param ops
 * @param para
 * @param workspace
 */
//获取P向量组
//V表示用于计算的基底向量组(GCGE_Vec)
void GCGE_ComputeP(GCGE_DOUBLE *subspace_evec, void **V, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
    //要给的形式参数
    GCGE_DOUBLE *coef    = subspace_evec;
    GCGE_INT    ldc      = workspace->dim_xpw;
    GCGE_INT    *unlock  = workspace->unlock;
    GCGE_INT    p_start  = workspace->dim_x;
    GCGE_INT    p_ncols  = workspace->unconv_bs;
    GCGE_INT    xx_ncols = workspace->dim_x;
    GCGE_INT    xx_nrows = workspace->last_dim_x;
    void        **V_tmp  = workspace->V_tmp;
    GCGE_INT    mv_s[2]  = { 0,0 };
    GCGE_INT    mv_e[2]  = { 0,0 };

    //函数中用到的参数定义
    GCGE_INT    current_p = 0;
    GCGE_INT    current_x = 0;
    GCGE_INT    p_end     = 0;
    GCGE_DOUBLE *px = coef + p_start * ldc;
    GCGE_DOUBLE *pp = px   + xx_nrows;
    GCGE_DOUBLE *xp = coef + xx_nrows;
    /*
    printf("221 subspace_evec: \n");
    printf("%f\t%f\t%f\n", coef[0], coef[3], coef[6]);
    printf("%f\t%f\t%f\n", coef[1], coef[4], coef[7]);
    printf("%f\t%f\t%f\n", coef[2], coef[5], coef[8]);
    */

    //初始化为0并进行复制
    memset(px, 0.0, p_ncols*ldc*sizeof(GCGE_DOUBLE));
    /*
    printf("228 subspace_evec: \n");
    printf("%f\t%f\t%f\n", coef[0], coef[3], coef[6]);
    printf("%f\t%f\t%f\n", coef[1], coef[4], coef[7]);
    printf("%f\t%f\t%f\n", coef[2], coef[5], coef[8]);
    */
    for(current_p=0; current_p<p_ncols; current_p++)
    {
        current_x = unlock[current_p];
        memcpy(pp+current_p*ldc, xp+current_x*ldc, (ldc-xx_nrows)*sizeof(GCGE_DOUBLE));
    }

    /*
    printf("subspace_evec: \n");
    printf("%f\t%f\t%f\n", coef[0], coef[3], coef[6]);
    printf("%f\t%f\t%f\n", coef[1], coef[4], coef[7]);
    printf("%f\t%f\t%f\n", coef[2], coef[5], coef[8]);
    */
    //子空间正交化
    p_end = p_start + p_ncols;
    //printf("p_start:%d, p_end: %d\n", p_start, p_end);
    GCGE_OrthogonalSubspace(coef, ldc, p_start, &p_end, NULL, -1, para->orth_para);
    //printf("p_start:%d, p_end: %d\n", p_start, p_end);
    p_ncols = p_end - p_start;

    //V线性组合得到P向量，存放在V_tmp中
    mv_s[0] = 0;
    mv_s[1] = 0;
    mv_e[0] = ldc;
    mv_e[1] = p_ncols;
    ops->MultiVecLinearComb(V, V_tmp, mv_s, mv_e, px, ldc, NULL, -1, ops);

    //把V_tmp与V中P的部分进行指针交换
    mv_s[0] = p_start;
    mv_s[1] = 0;
    mv_e[0] = p_end;
    mv_e[1] = p_ncols;
    ops->MultiVecSwap(V, V_tmp, mv_s, mv_e, ops);

    //更新workspace中的dim_xp
    workspace->dim_xp = p_end;

}

