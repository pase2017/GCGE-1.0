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

/* brief: 计算Rayleigh-Ritz过程中的子空间矩阵V^T*A*V中的大规模操作
 *        subspace_matrix = [X,P,W]^T * A * [W] = XAW
 *                                                PAW
 *                                                WAW
 *
 *   需要 submat, ldm, W_start, W_end, V_tmp(临时存储A*W)
 *
 */
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
    //mv_s,mv_e表示每个多向量操作函数中，两个向量组中操作向量的起始和终点位置
    
    //比如下面表示对V的mv_s[0]到mv_e[0]-1位置的向量进行操作
    //    对workspace的mv_s[1]到mv_e[1]-1位置的向量进行操作
    
    //计算 workspace(0:w_length) = A * W(start[1]:end[1])
    //这里workspace 为大规模向量组临时空间
    ops->MatDotMultiVec(A, V, workspace, mv_s, mv_e, ops);

    mv_s[0] = 0;
    mv_e[0] = xp_length;
    mv_s[1] = 0;
    mv_e[1] = w_length;

    //计算 subspace_mat = V * workspace
    //这里计算 X(AW) 部分
    //         P(AW)
    //需要 subspace_mat 的 leading dimension, 这里为end_V
    ops->MultiVecInnerProd(V, workspace, subspace_mat, "nonsym", mv_s, mv_e, end_V, ops);

    //double *vTAw = subspace_mat;

    mv_s[0] = start_W;
    mv_e[0] = end_V;
    mv_s[1] = 0;
    mv_e[1] = w_length;

    //计算 subspace_mat = W * workspace, 且已知 subspace_mat 是对称矩阵
    //这里计算 W(AW) 部分
    //需要 subspace_mat 的 leading dimension, 这里为end_V
    //     以及子空间矩阵存放的起始位置 subspace_mat+xp_length
    ops->MultiVecInnerProd(V, workspace, subspace_mat + xp_length, "sym", mv_s, mv_e, end_V, ops);

}


//计算子空间矩阵subspace_matrix=V^TAV
//A是计算子空间矩阵时用到的矩阵Ａ(GCGE_MATRIX), V是用到的向量组(GCGE_Vec)
//void GCGE_GetSubspaceMatrix(void *A, void **V, GCGE_DOUBLE *subspace_matrix, 
/* brief: 计算Rayleigh-Ritz过程中的子空间矩阵V^T*A*V
 *        subspace_matrix = [X,P,W]^T * A * [X,P,W] = XAX  XAP  XAW
 *                                                    PAX  PAP  PAW
 *                                                    WAX  WAP  WAW
 *   1. 小规模操作,需要 coef, ldc, AA_last, p_end, 
 *                      ldm(新矩阵subspace_matrix的规模leading dimension)
 *
 *        XAX  XAP  = XX  PX ^T * AA_last *  XX  PX
 *        PAX  PAP    XP  PP                 XP  PP
 *                    XW  PW                 XW  PW
 *
 *   2. 大规模操作，计算 XAW , 需要 submat, ldm, W_start, W_length,
 *                       PAW        V_tmp(临时存储A*W)
 *                       WAW
 *
 */
/* TODO need change input */
void GCGE_ComputeSubspaceMatrix(void *A, void **V, 
        GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
     //ldm: laplack的风格名字：表示subspace_matrix的维数
    GCGE_INT    i, dim_xp = workspace->dim_xp, ldm = workspace->dim_xpw;
    GCGE_DOUBLE *subspace_dtmp   = workspace->subspace_dtmp; 
    /*subspace_dtmp表示下式左端矩阵的起始位置
     *  XAX  XAP  = XX  PX ^T * AA_last *  XX  PX
     *  PAX  PAP    XP  PP                 XP  PP
     *              XW  PW                 XW  PW 
     */
    GCGE_DOUBLE *subspace_matrix = workspace->subspace_matrix;
    GCGE_DOUBLE *subspace_evec   = workspace->subspace_evec;
    /* TODO */
    if(dim_xp != 0)
    {
    	    //lde表示AA_last的行数(leading dimension of subspace_evec)
    	    //AA_last: 上一次迭代的AA矩阵的行数
        GCGE_INT    lde = workspace->last_dim_xpw; 
        GCGE_DOUBLE alpha = 1.0, beta = 0.0;
        // 这里只是取了一个临时空间用于存储下式
        // AP表示 AA_last *  XX  PX 的起始位置, 
        //                   XP  PP
        //                   XW  PW
        // AA_last此时存储在subspace_matrix中
        // dsymm 计算AP,因为AA_last是对称矩阵，所以使用dsymm计算对称矩阵乘矩阵
        // 需要的参数有 AA_last 的行数 lde，XX PX 的总列数 dim_xp
        //DenseMatDotDenseMat的三个矩阵中是否有可覆盖的可能（为了减少内存开销）
        ops->DenseSymMatDotDenseMat("L", "U", &lde, &dim_xp, &alpha, subspace_matrix, 
                &lde, workspace->subspace_evec, &lde, &beta, subspace_dtmp, &lde);
        //dgemm 计算下式
        //  XAX  XAP  = XX  PX ^T * AP
        //  PAX  PAP    XP  PP    
        //              XW  PW    
        // 需要的参数有 左矩阵的行数 lde, 右矩阵的列数 dim_xp, 
        // 左矩阵的列数(右行) lde (这里需要添加几个参数命名 TODO)
        memset(subspace_matrix, 0.0, ldm*ldm*sizeof(GCGE_DOUBLE));
        ops->DenseMatDotDenseMat("T", "N", &dim_xp, &dim_xp, &lde, &alpha, 
                workspace->subspace_evec, &lde, subspace_dtmp, 
                &lde, &beta, subspace_matrix, &ldm);

    }
    GCGE_ComputeSubspaceMatrixVTAW(V, A, dim_xp, ldm, subspace_matrix+dim_xp*ldm, 
	  ops, workspace->evec);
}

/* brief: 计算subspace_matrix的特征对
 *    
 *    需要的参数：subspace_matrix的规模: ldm
 *                要计算特征值的范围：il,iu(计算第il到第iu个特征对)
 *
 *    这里 ops->DenseMatEigenSolver 的形式参数给了 dsyevx和dsyevr 用到的
 *    所有参数，工作空间给了两个函数需要的最大空间，
 *    默认选择dsyevx，可以在使用时设置使用哪一种求解方式
 */
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
    //以下N代表ldm
    //work给足够大的空间: max(26*N, 1+6*N+2*N**2)
    //lwork: dsyevr(26*N),dsyevx(8*N),dsyev(3*N),dsyevd(1+6*N+2*N**2)
    //workspace_dtmp空间需要；N*N+lwork
    //iwork足够大的空间：10*N
    //liwork: dysevr(10*N),dsyevx(5*N,没有这个参数),
    //        dsyev(没有这个参数，且不需要iwork),dsyevd(3+5*N)
    //isuppz: dsyevr(2*N的空间)
    //ifail: dsyevx(N的空间）
    ops->DenseMatEigenSolver("V", "I", "U", &nrows, temp_matrix, &ncols, 
            &vl, &vu, &il, &iu, &abstol, 
            &m, eval, subspace_evec, &ldm, 
            isuppz, dwork_space, &lwork,
            subspace_itmp, &liwork, ifail, &info);
    if(para->use_mpi_bcast)
    {
#if GCGE_USE_MPI
        MPI_Bcast(subspace_evec, m*ldm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    }
    /*
    GCGE_INT i, j, k;
    for(i=il-1; i<iu; i++)
    {
        if(GCGE_ModuleMaxDouble(subspace_evec+i*ldm, ldm) < 0.0)
        {
            for(j=0; j<ldm; j++)
            {
                subspace_evec[i*ldm+j] *= -1.0;
            }
        }
    }
    */
    /* 为了保证不同进程间特征向量的方向一致，
     * 取每个特征向量的最后一个非零元素为正数(若不是，则这个特征向量乘以-1)
     * 一般做法是把每一行的绝对值最大值的位置变成正数
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
            else if(subspace_evec[i*ldm+j] > -negative_zero_tol)
            { 
               //正数足够大之后可以结束本向量的搜索
               break;
            }
        }
    }
     */
    //用lapack_syev计算得到的特征值是按从小到大排列,
    //如果用A内积，需要把特征值取倒数后再按从小到大排列
    //由于后面需要用到的只有前dim_x个特征值，所以只把前dim_x个拿到前面来
    if(strcmp(para->orth_type, "A") == 0)
    {
        GCGE_SortEigenpairs(eval, subspace_evec, iu, ldm, workspace->subspace_dtmp);
    }
}

//用A内积的话，特征值从小到大，就是原问题特征值倒数的从小到大，
//所以顺序反向，同时特征值取倒数
void GCGE_SortEigenpairs(GCGE_DOUBLE *eval, GCGE_DOUBLE *evec, GCGE_INT nev, 
					GCGE_INT ldv, GCGE_DOUBLE *work)
{
    GCGE_INT head = 0, tail = ldv-1;
    GCGE_DOUBLE tmp;
    for( head=0; head<nev; head++ )
    {
        tail = ldv-1-head;
        if(head < tail)
        {
            memcpy(work, evec+head*ldv, ldv*sizeof(GCGE_DOUBLE));
            memcpy(evec+head*ldv, evec+tail*ldv, ldv*sizeof(GCGE_DOUBLE));
            memcpy(evec+tail*ldv, work, ldv*sizeof(GCGE_DOUBLE));
            tmp = eval[head];
            eval[head] = 1.0/eval[tail];
            eval[tail] = 1.0/tmp;
        }
        else
        {
            break;
        }
    }
}

/*
GCGE_DOUBLE GCGE_ModuleMaxDouble(GCGE_DOUBLE *a, GCGE_INT n)
{
    GCGE_INT i = 0;
    GCGE_DOUBLE max = a[0];
    for(i=1; i<n; i++)
    {
        if(fabs(max+ 1e-16) < fabs(a[i]))
            max = a[i];
    }
    return max;
}

*/
