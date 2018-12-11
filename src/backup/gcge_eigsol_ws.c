/*
 * =====================================================================================
 *
 *       Filename:  gcge_eigsol_ws.c
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

#include "gcge_eigsol_ws.h"

//给GCGE_EIGSOL_WS结构中各个部分分配空间
void GCGE_EIGSOL_WS_Create(GCGE_EIGSOL_WS **eigsol_ws, GCGE_EIGSOL_PARA *para, 
		GCGE_OPS *ops, void *A)
{
    (*eigsol_ws) = (GCGE_EIGSOL_WS *)malloc(sizeof(GCGE_EIGSOL_WS));

	//RitzVec临时存储Ｘ, X的个数上限为nev加重数
	GCGE_INT size_RitzVec = (GCGE_INT)(para->nev * (1.0+para->eval_multi_max));
	//V存储[X,P,W],P,W的个数上限为block_size
	GCGE_INT size_V       = size_RitzVec + 2 * para->block_size;
	//V_tmp在正交化中(块:block_size,GS:1), RR中(VTAW:block_size),
	//ComputeP中(暂时存储Ｐ:block_size),ComputeW(rhs:1),CG(3)
	GCGE_INT size_V_tmp   = (para->block_size > 4)?(para->block_size):4;

	//子空间矩阵存储VTAV,大小上限与Ｖ的空间上限有关
	GCGE_INT size_submat  = size_V * size_V;
	//subeval需要算X的上限个数
	GCGE_INT size_subeval = size_RitzVec;
	//subevec需要算Ｘ的上限个数,进行Ｐ的计算时需要再加Ｐ的个数列
	GCGE_INT size_subevec = size_V * (size_V + para->block_size);
	//dsyevr: liwork(10*N),isuppz:(2*N); dsyevx: liwork(5*N),ifail(N)
	GCGE_INT size_subitmp = 12 * size_V;
    //lwork: dsyevr(26*N,实际测试33*N),dsyevx(8*N,实际测试更多),dsyev(3*N),dsyevd(1+6*N+2*N**2)
	//临时存储submat
	GCGE_INT size_subdtmp = 33 * size_V + 3 * size_V * size_V;

	(*eigsol_ws)->space_size_RitzVec = size_RitzVec;
	(*eigsol_ws)->space_size_V       = size_V;
	(*eigsol_ws)->space_size_V_tmp   = size_V_tmp;
    
	(*eigsol_ws)->unlock = (GCGE_INT*)calloc(para->nev, sizeof(GCGE_DOUBLE));
	(*eigsol_ws)->subspace_matrix = (GCGE_DOUBLE*)calloc(size_submat,  sizeof(GCGE_DOUBLE));
	(*eigsol_ws)->subspace_evec   = (GCGE_DOUBLE*)calloc(size_subevec, sizeof(GCGE_DOUBLE));
	(*eigsol_ws)->subspace_eval   = (GCGE_DOUBLE*)calloc(size_subeval, sizeof(GCGE_DOUBLE));
	(*eigsol_ws)->subspace_dtmp   = (GCGE_DOUBLE*)calloc(size_subdtmp, sizeof(GCGE_DOUBLE));
	(*eigsol_ws)->subspace_itmp   = (GCGE_INT*)calloc(size_subitmp, sizeof(GCGE_DOUBLE));

    ops->BuildMultiVecByMat(A, &((*eigsol_ws)->V),       size_V,       ops);
    ops->BuildMultiVecByMat(A, &((*eigsol_ws)->V_tmp),   size_V_tmp,   ops);
    ops->BuildMultiVecByMat(A, &((*eigsol_ws)->RitzVec), size_RitzVec, ops);
}

//释放GCGE_Para结构中分配的空间
void GCGE_EIGSOL_WS_Free(GCGE_EIGSOL_WS **eigsol_ws, GCGE_OPS *ops)
{
    if((*eigsol_ws)->unlock)
    {
        free((*eigsol_ws)->unlock);
        (*eigsol_ws)->unlock = NULL;
    }
    if((*eigsol_ws)->subspace_matrix)
    {
        free((*eigsol_ws)->subspace_matrix);
        (*eigsol_ws)->subspace_matrix = NULL;
    }
    if((*eigsol_ws)->subspace_evec)
    {
        free((*eigsol_ws)->subspace_evec);
        (*eigsol_ws)->subspace_evec = NULL;
    }
    if((*eigsol_ws)->subspace_eval)
    {
        free((*eigsol_ws)->subspace_eval);
        (*eigsol_ws)->subspace_eval = NULL;
    }
    if((*eigsol_ws)->subspace_dtmp)
    {
        free((*eigsol_ws)->subspace_dtmp);
        (*eigsol_ws)->subspace_dtmp = NULL;
    }
    if((*eigsol_ws)->subspace_itmp)
    {
        free((*eigsol_ws)->subspace_itmp);
        (*eigsol_ws)->subspace_itmp = NULL;
    }
    if((*eigsol_ws)->V)
    {
        ops->MultiVecDestroy(&((*eigsol_ws)->V), (*eigsol_ws)->space_size_V, ops);
    }
    if((*eigsol_ws)->V)
    {
        ops->MultiVecDestroy(&((*eigsol_ws)->V_tmp), (*eigsol_ws)->space_size_V_tmp, ops);
    }
    if((*eigsol_ws)->V)
    {
        ops->MultiVecDestroy(&((*eigsol_ws)->RitzVec), (*eigsol_ws)->space_size_RitzVec, ops);
    }

    free((*eigsol_ws)); (*eigsol_ws) = NULL;
}
