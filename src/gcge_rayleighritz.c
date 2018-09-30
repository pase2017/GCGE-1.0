/*
 * =====================================================================================
 *
 *       Filename:  gcge_rayleighritz.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gcge_rayleighritz.h"

//计算subspace_matrix = V^T*A*W
void GCGE_ComputeSubspaceMatrixVTAW(void **V, void *A, GCGE_INT start_W, GCGE_INT end_V, GCGE_DOUBLE *subspace_mat, GCGE_OPS *ops, void **workspace)
{
    //start[0],end[0]表示V中[X,P]的起始与终点位置
    GCGE_INT w_length = end_V - start_W;
    GCGE_INT xp_length = end_V - w_length;

    GCGE_INT mv_s[2];
    GCGE_INT mv_e[2];

    mv_s[0] = start_W;
    mv_e[0] = end_V;
    mv_s[1] = 0;
    mv_e[1] = w_length;
    
    printf("mv_s: %d, %d, mv_e: %d, %d\n", mv_s[0], mv_s[1], mv_e[0], mv_e[1]);
    //计算 workspace(0:w_length) = A * W(start[1]:end[1])
    ops->MatDotMultiVec(A, V, workspace, mv_s, mv_e, ops);

    mv_s[0] = 0;
    mv_e[0] = xp_length;
    mv_s[1] = 0;
    mv_e[1] = w_length;

    //计算 subspace_mat = V * workspace
    ops->MultiVecInnerProd(V, workspace, subspace_mat, "nonsym", mv_s, mv_e, end_V, ops);

    double *vTAw = subspace_mat;

//    printf ( "VTAW\n" );


    mv_s[0] = start_W;
    mv_e[0] = end_V;
    mv_s[1] = 0;
    mv_e[1] = w_length;

    //计算 subspace_mat = W * workspace, 且已知 subspace_mat 是对称矩阵
    ops->MultiVecInnerProd(V, workspace, subspace_mat + xp_length, "sym", mv_s, mv_e, end_V, ops);

//    printf ( "%f\t%f\n",  vTAw[0],  vTAw[3] );
//    printf ( "%f\t%f\n",  vTAw[1],  vTAw[4] );
//    printf ( "%f\t%f\n",  vTAw[2],  vTAw[5] );

}


//计算子空间矩阵subspace_matrix=V^TAV
//A是计算子空间矩阵时用到的矩阵Ａ(GCGE_MATRIX), V是用到的向量组(GCGE_Vec)
//void GCGE_GetSubspaceMatrix(void *A, void **V, GCGE_DOUBLE *subspace_matrix, 
/* TODO need change input */
void GCGE_ComputeSubspaceMatrix(void *A, void **V, 
        GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
    GCGE_INT    i, dim_xp = workspace->dim_xp, ldm = workspace->dim_xpw;
    GCGE_DOUBLE *subspace_dtmp   = workspace->subspace_dtmp;
    GCGE_DOUBLE *subspace_matrix = workspace->subspace_matrix;
    GCGE_DOUBLE *subspace_evec   = workspace->subspace_evec;
    /* TODO */
#if 0
    t1 = GCGE_GetTime();
#endif
    if(dim_xp != 0)
    {
        GCGE_DOUBLE alpha = 1.0, beta = 0.0,
                    *AP = subspace_dtmp+dim_xp*dim_xp;
        GCGE_INT    lde = workspace->last_dim_xpw;
        //dsymm 计算AP = submat * subevec
        ops->DenseSymMatDotDenseMat("L", "U", &lde, &dim_xp, &alpha, subspace_matrix, 
                &lde, workspace->subspace_evec, &lde, &beta, AP, &lde);
        //dgemm 计算submat = subevec^T * AP
        ops->DenseMatDotDenseMat("T", "N", &dim_xp, &dim_xp, &lde, &alpha, 
                workspace->subspace_evec, &lde, AP, 
                &lde, &beta, subspace_dtmp, &dim_xp);
        memset(subspace_matrix, 0.0, ldm*ldm*sizeof(GCGE_DOUBLE));
        for( i=0; i<dim_xp; i++ )
        {
            memcpy(subspace_matrix+i*ldm, subspace_dtmp+i*dim_xp, 
                    dim_xp*sizeof(GCGE_DOUBLE));
        }
    }
#if 0
    GCGE_DOUBLE t3 = GCGE_GetTime();
#endif
    printf("dim_xp: %d\n", dim_xp);
    GCGE_ComputeSubspaceMatrixVTAW(V, A, dim_xp, ldm, subspace_matrix+dim_xp*ldm, 
	  ops, workspace->V_tmp);
#if 0
    t2 = GCGE_GetTime();
    sub1 += t3-t1;
    sub2 += t2-t3;
    //统计RayleighRitz问题求解时间
    para->stat_para->PartTimeTotal->RayleighRitzMat_Time += t2-t1;
    para->stat_para->PartTimeInOneIte->RayleighRitzMat_Time = t2-t1;
#endif
}

//用LAPACKE_dsyev求解子空间矩阵的特征对
void GCGE_ComputeSubspaceEigenpairs(GCGE_DOUBLE *subspace_matrix, 
        GCGE_DOUBLE *eval, GCGE_DOUBLE *subspace_evec, 
        GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
    /* TODO */
#if 0
    GCGE_DOUBLE    t1 = GCGE_GetTime();
#endif
    GCGE_INT    ldm = workspace->dim_xpw;
    GCGE_INT    info;
    GCGE_DOUBLE *temp_matrix = workspace->subspace_dtmp,
                *dwork_space = temp_matrix + ldm*ldm,
                vl = 0.0, vu = 0.0, abstol = 0.0;
    //abstol = 2*dlamch_("S");
    GCGE_INT    max_dim_x = workspace->max_dim_x;
    GCGE_INT    il = 1, 
                iu = (max_dim_x < ldm) ? max_dim_x : ldm,
                m = iu-il+1;
    GCGE_INT    lwork = 8*ldm;
    GCGE_INT    liwork = 5*ldm;
    GCGE_INT    *subspace_itmp = workspace->subspace_itmp;
    //subspace_itmp[liwork]在dsyevr时存储liwork
    GCGE_INT    *isuppz = subspace_itmp + 10*ldm; 
    GCGE_INT    *ifail  = subspace_itmp + 5*ldm; 
    memcpy(temp_matrix, subspace_matrix, ldm*ldm*sizeof(GCGE_DOUBLE));
    
    GCGE_INT    nrows = ldm, ncols = ldm;
    //printf("ldm: %d\n", ldm);
    //printf("dwork_sapce: %p, %p\n", dwork_space, dwork_space+lwork);
    //以下N代表ldm
    //work给足够大的空间: max(26*N, 1+6*N+2*N**2)
    //lwork: dsyevr(26*N),dsyevx(8*N),dsyev(3*N),dsyevd(1+6*N+2*N**2)
    //workspace_dtmp空间需要；N*N+lwork
    //iwork足够大的空间：10*N
    //liwork: dysevr(10*N),dsyevx(5*N,没有这个参数),
    //        dsyev(没有这个参数，且不需要iwork),dsyevd(3+5*N)
    //isuppz: dsyevr(2*N的空间)
    //ifail: dsyevx(N的空间）
    /*
    printf("num_iter: %d\n", para->num_iter);
    if(para->num_iter > 1)
    {
        lwork = -1;
    }
    */
    ops->DenseMatEigenSolver("V", "I", "U", &nrows, temp_matrix, &ncols, 
            &vl, &vu, &il, &iu, &abstol, 
            &m, eval, subspace_evec, &ldm, 
            isuppz, dwork_space, &lwork,
            subspace_itmp, &liwork, ifail, &info);
    //printf("ldm: %d, lwork: %d, work[0]:  %f\n", ldm, lwork, dwork_space[0]);
    //使所有会用到的特征向量的最后一个非0分量保持为正,避免出现不同进程特征向量方向不一致的问题
    GCGE_INT i, j, k;
    GCGE_DOUBLE negative_zero_tol = -1e-10;
    for(i=0; i< iu; i++)
    {
        for(j=ldm-1; j>0; j--)
        {
            if(subspace_evec[i*ldm+j] < negative_zero_tol)
            {
               for(k=0; k<ldm; k++)
               {
                  subspace_evec[i*ldm+k] *= -1.0;
               }
               break;
            }
            else if(subspace_evec[i*ldm+j] > negative_zero_tol)
            {
               break;
            }
        }
    }
    //用lapack_syev计算得到的特征值是按从小到大排列,
    //如果用A内积，需要把特征值取倒数后再按从小到大排列
    //由于后面需要用到的只有前dim_x个特征值，所以只把前dim_x个拿到前面来
    if(strcmp(para->orth_type, "A") == 0)
    {
        GCGE_SortEigenpairs(eval, subspace_evec, iu, ldm, workspace->subspace_dtmp);
    }
#if 0
    GCGE_DOUBLE t2 = GCGE_GetTime();
    //统计RayleighRitz问题求解时间
    para->stat_para->PartTimeTotal->RayleighRitzEigen_Time += t2-t1;
    para->stat_para->PartTimeInOneIte->RayleighRitzEigen_Time = t2-t1;
#endif
}

//用A内积的话，特征值从小到大，就是原问题特征值倒数的从小到大，
//所以顺序反向，同时特征值取倒数
void GCGE_SortEigenpairs(GCGE_DOUBLE *eval, GCGE_DOUBLE *evec, GCGE_INT nev, GCGE_INT ldv, GCGE_DOUBLE *work)
{
    GCGE_INT head = 0, tail = ldv-1;
    for( head=0; head<nev; head++ )
    {
        tail = ldv-1-head;
        if(head < tail)
        {
            memcpy(work, evec+head*ldv, ldv*sizeof(GCGE_DOUBLE));
            memcpy(evec+head*ldv, evec+tail*ldv, ldv*sizeof(GCGE_DOUBLE));
            memcpy(evec+tail*ldv, work, ldv*sizeof(GCGE_DOUBLE));
            work[0] = eval[head];
            eval[head] = 1.0/eval[tail];
            eval[tail] = 1.0/work[0];
        }
        else
        {
            break;
        }
    }
}


