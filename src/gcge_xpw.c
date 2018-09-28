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


//计算RitzVec = V*subspace_evec
//V是用来计算Ritz向量的基底向量组(GCGE_Vec),RitzVec是计算得到的Ritz向量(GCGE_Vec)
void GCGE_ComputeRitzVectors(void **V, GCGE_DOUBLE *subspace_evec, void **RitzVec, 
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

//用未收敛的特征向量作为X求解线性方程组 A W = \lambda B X
//A,B表示矩阵A,B(GCGE_MATRIX),V表示向量组(GCGE_Vec)
void GCGE_ComputeW(void *A, void *B, void **V, GCGE_DOUBLE *eval, 
        GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
    /* TODO */
#if 0
    GCG_Double    t1 = GCG_GetTime();
#endif
    GCGE_INT i, j, k, start = workspace->dim_xp,
            *unlock = workspace->unlock;
    void    **V_tmp = workspace->V_tmp;
    for( i=0; i < workspace->unconv_bs; i++ )
    {
        j = unlock[i];
        //初值V[start+i]=V[i]
        //Vtmp[0]=\lambda*B*V[i]作为右端项
        //计算V[start+i]=A^(-1)BV[i]
        //调用CG迭代来计算线性方程组
        if(para->if_lobgcg == 0)
        {
            //使用GCG
            ops->VecAxpby(1.0, V[j], 0.0, V[start+i]);
            if(B == NULL)
            {
                ops->VecAxpby(eval[j], V[j], 0.0, V_tmp[0]);
            }
            else
            {
                ops->MatDotVec(B, V[j], V_tmp[0]);
                ops->VecAxpby(0.0, V_tmp[0], eval[j], V_tmp[0]);
            }
        }
        else
        {
            //使用LOBGCG
			ops->MatDotVec(A, V[j], V_tmp[0]);
            if(B == NULL)
            {
                ops->VecAxpby(-eval[j], V[j], 1.0, V_tmp[0]);
            }
            else
            {
			    ops->MatDotVec(B, V[j], V_tmp[1]);
                ops->VecAxpby(-eval[j], V_tmp[1], 1.0, V_tmp[0]);
            }
            ops->VecAxpby(0.0, V[start+i], 0.0, V[start+i]);
        }
        if(para->if_use_cg == 1)
        {
            GCGE_CG(A, V_tmp[0], V[start+i], ops, para, V_tmp);
        }
        else
        {
            ops->LinearSolver(A, V_tmp[0], V[start+i], ops);
        }
    }
    workspace->last_dim_xpw = workspace->dim_xpw;
    workspace->dim_xpw = workspace->dim_xp + workspace->unconv_bs;
#if 0
    GCG_Double t2 = GCG_GetTime();
    para->stat_para->PartTimeTotal->GetW_Time += t2-t1;
    para->stat_para->PartTimeInOneIte->GetW_Time = t2-t1;
#endif
}

//获取P向量组
//V表示用于计算的基底向量组(GCGE_Vec)
void GCGE_ComputeP(GCGE_DOUBLE *subspace_evec, void **V, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
    /* TODO */
#if 0
    GCGE_DOUBLE    t1 = GCGE_GetTime();
#endif
    GCGE_INT i, j, n_vec, dim_x = workspace->dim_x,
             last_dim_x = workspace->last_dim_x,
             dim_xpw = workspace->dim_xpw,
             *unlock = workspace->unlock;
    void     **V_tmp = workspace->V_tmp;
#if 0
    solver->para->stat_para->ite_time = t1;
#endif
    //小规模正交化，构造d,构造subspace_matrix_sub用于下次计算
#if !ORTHPEFF
    for( i=0; i< workspace->unconv_bs; i++ )
    {
        memset(subspace_evec+(dim_x+i)*dim_xpw, 0.0, dim_xpw*sizeof(GCGE_DOUBLE));
        memcpy(subspace_evec+(dim_x+i)*dim_xpw+last_dim_x, 
               subspace_evec+unlock[i]*dim_xpw+last_dim_x, 
               (dim_xpw-last_dim_x)*sizeof(GCGE_DOUBLE));
    }
#endif
    //GCGE_DOUBLE test_time2 = GCGE_GetTime();
    //更新dim_xp,dim_xp表示[X,P]的向量个数
    workspace->dim_xp = workspace->dim_x+workspace->unconv_bs;
    GCGE_OrthogonalSubspace(subspace_evec, dim_x, &(workspace->dim_xp), dim_xpw, NULL, 0, para->orth_para, workspace);
#if 0
    GCGE_DOUBLE t2 = GCGE_GetTime();
    solver->para->stat_para->PartTimeTotal->GetPOrth_Time += t2-t1;
    solver->para->stat_para->PartTimeInOneIte->GetPOrth_Time = t2-t1;
#endif
    //用基向量组V计算P所对应的长向量，存在V_tmp中
    n_vec = workspace->dim_xp - dim_x;
    GCGE_INT start[2] = {0, 0};
    GCGE_INT end[2]   = {dim_xpw, n_vec};

    if(ops->MultiVecLinearComb)
    {
        ops->MultiVecLinearComb(V, V_tmp, start, end, subspace_evec+dim_x*dim_xpw, dim_xpw, NULL, 0, ops);

        start[0] = dim_x;
        end[0]   = workspace->dim_xp;
        start[1] = 0;
        end[1]   = n_vec;

        ops->MultiVecSwap(V, V_tmp, start, end, ops);
    }
    //交换V_tmp与V[dim_x:dim_xp]的向量指针，使Ｐ存到V的dim_x到dim_xp列中
#if 0
    GCGE_DOUBLE t3 = GCGE_GetTime();
    solver->para->stat_para->PartTimeTotal->GetPLinearAlg_Time += t3-t2;
    solver->para->stat_para->PartTimeInOneIte->GetPLinearAlg_Time = t3-t2;
#endif
}

