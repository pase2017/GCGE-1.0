/*
 * =====================================================================================
 *
 *       Filename:  gcge_workspace.c
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

#include "gcge_workspace.h"

//给GCGE_Para结构中各个部分分配空间
void GCGE_WORKSPACE_Create(GCGE_WORKSPACE **workspace)
{
    (*workspace) = (GCGE_WORKSPACE*)malloc(sizeof(GCGE_WORKSPACE));
    //对GCGE算法中用到的一些参数进行初始化
    (*workspace)->max_dim_x = 0;
    (*workspace)->dim_x     = 0;
    (*workspace)->dim_xp    = 0;
    (*workspace)->dim_xpw   = 0;
    (*workspace)->conv_bs   = 0;
    (*workspace)->unconv_bs = 0;

    //近似特征值
    (*workspace)->eval          = NULL;
    //小规模的临时工作空间
    (*workspace)->subspace_dtmp = NULL;
    (*workspace)->subspace_itmp = NULL;
    //存储子空间矩阵的特征向量
    (*workspace)->subspace_evec = NULL;
    //用于存储子空间矩阵
    (*workspace)->subspace_matrix = NULL;

    //存储前nev个特征对中未收敛的特征对编号
    (*workspace)->unlock = NULL;
    //正交化时用到的临时GCGE_INT*型空间,用于临时存储非0的编号
    (*workspace)->orth_ind = NULL;
    (*workspace)->V = NULL;
    (*workspace)->V_tmp = NULL;
    (*workspace)->RitzVec = NULL;
    (*workspace)->CG_p = NULL;
    (*workspace)->CG_r = NULL;
    
    //GCGE_STATISTIC_Para用于统计各部分时间
}
void GCGE_WORKSPACE_Setup(GCGE_WORKSPACE *workspace, GCGE_PARA *para, GCGE_OPS *ops, void *A)
{
    GCGE_INT       i, nev    = para->nev,
                   max_dim_x = (nev*1.25 > nev+2) ? (nev*1.25) : (nev+2);

    //对GCGE算法中用到的一些参数进行初始化
    workspace->max_dim_x = max_dim_x;
    workspace->dim_x     = nev;
    workspace->dim_xp    = 0;
    workspace->dim_xpw   = nev;
    workspace->conv_bs   = 0;
    workspace->unconv_bs = (nev < para->block_size) ? nev : para->block_size;

    GCGE_INT max_dim_xpw = max_dim_x + 2 * para->block_size;
    //近似特征值
    workspace->eval          = (GCGE_DOUBLE*)calloc(max_dim_xpw, sizeof(GCGE_DOUBLE));
    //小规模的临时工作空间
    //GCGE_INT lwork1 = 26*max_dim_xpw;
    //GCGE_INT lwork2 = 1+6*max_dim_xpw+2*max_dim_xpw*max_dim_xpw;
    workspace->subspace_dtmp = (GCGE_DOUBLE*)calloc(max_dim_xpw*max_dim_xpw+300*max_dim_xpw, sizeof(GCGE_INT));
            //+(lwork1>lwork2)?lwork1:lwork2, sizeof(GCGE_DOUBLE));
    workspace->subspace_itmp = (GCGE_INT*)calloc(12*max_dim_xpw, sizeof(GCGE_INT));
    //存储子空间矩阵的特征向量
    workspace->subspace_evec = (GCGE_DOUBLE*)calloc(max_dim_xpw*max_dim_xpw, sizeof(GCGE_DOUBLE));
    //用于存储子空间矩阵
    workspace->subspace_matrix = (GCGE_DOUBLE*)calloc(max_dim_xpw*max_dim_xpw, sizeof(GCGE_DOUBLE));

    //存储前nev个特征对中未收敛的特征对编号
    workspace->unlock = (GCGE_INT*)calloc(max_dim_x, sizeof(GCGE_INT));
    //for(i=0; i<nev; i++)
    for(i=0; i<max_dim_x; i++)
    {
        workspace->unlock[i] = i;
    }
    //正交化时用到的临时GCGE_INT*型空间,用于临时存储非0的编号
    workspace->orth_ind = (GCGE_INT*)calloc(max_dim_xpw, sizeof(GCGE_INT));

    //V,V_tmp,RitzVec是向量工作空间
    ops->BuildMultiVecByMat(A, &(workspace->V), max_dim_xpw, ops);
    ops->BuildMultiVecByMat(A, &(workspace->V_tmp), (max_dim_xpw>4)?max_dim_xpw:4, ops);
    ops->BuildMultiVecByMat(A, &(workspace->RitzVec), max_dim_x, ops);
    ops->BuildMultiVecByMat(A, &(workspace->CG_p), para->block_size, ops);
    ops->BuildMultiVecByMat(A, &(workspace->CG_r), para->block_size, ops);

    //GCGE_STATISTIC_Para用于统计各部分时间
}
//释放GCGE_Para结构中分配的空间
void GCGE_WORKSPACE_Free(GCGE_WORKSPACE **workspace, GCGE_PARA *para, GCGE_OPS *ops)
{
    if((*workspace)->unlock)
    {
        free((*workspace)->unlock);
        (*workspace)->unlock = NULL;
    }
    if((*workspace)->subspace_matrix)
    {
        free((*workspace)->subspace_matrix);
        (*workspace)->subspace_matrix = NULL;
    }
    if((*workspace)->subspace_evec)
    {
        free((*workspace)->subspace_evec);
        (*workspace)->subspace_evec = NULL;
    }
    if((*workspace)->eval)
    {
        free((*workspace)->eval);
        (*workspace)->eval = NULL;
    }
    if((*workspace)->subspace_dtmp)
    {
        free((*workspace)->subspace_dtmp);
        (*workspace)->subspace_dtmp = NULL;
    }
    if((*workspace)->subspace_itmp)
    {
        free((*workspace)->subspace_itmp);
        (*workspace)->subspace_itmp = NULL;
    }
    if((*workspace)->orth_ind)
    {
        free((*workspace)->orth_ind);
        (*workspace)->orth_ind = NULL;
    }
    
    GCGE_INT max_dim_x = (*workspace)->max_dim_x;
    GCGE_INT max_dim_xpw = max_dim_x + 2 * para->block_size;
    if((*workspace)->V)
    {
        ops->FreeMultiVec(&((*workspace)->V), max_dim_xpw, ops);
        ops->FreeMultiVec(&((*workspace)->V_tmp), (max_dim_xpw>4)?max_dim_xpw:4, ops);
        ops->FreeMultiVec(&((*workspace)->RitzVec), max_dim_x, ops);
        ops->FreeMultiVec(&((*workspace)->CG_p), para->block_size, ops);
        ops->FreeMultiVec(&((*workspace)->CG_r), para->block_size, ops);
        if(ops->DenseMatCreate)
        {
            ops->DenseMatDestroy(&((*workspace)->dense_matrix));
        }
    }

    free((*workspace)); (*workspace) = NULL;
}
