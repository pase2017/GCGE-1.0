/*
 * =====================================================================================
 *
 *       Filename:  gcge_workspace.h
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
#ifndef _GCGE_WORKSPACE_H_
#define _GCGE_WORKSPACE_H_

#include "gcge_type.h"

#include "gcge_para.h"
#include "gcge_ops.h"

typedef struct GCGE_WORKSPACE_ {
    //向量工作空间
    void       **V; //V用于存储[X,P,W]基底
    void       **V_tmp; //V_tmp用作临时的向量存储空间
    /* TODO */
    void       **RitzVec; //RitzVec用于存储Ritz向量

    void       **CG_p; //CG迭代中的p向量
    void       **CG_r; //CG迭代中的r向量

    GCGE_INT    *unlock; //unlock用于存储未锁定的特征对编号(0-nev)
    GCGE_INT    max_dim_x; //X中最多的向量个数
    GCGE_INT    dim_x; //X中的向量个数
    GCGE_INT    dim_xp; //[X,P]中的向量个数
    GCGE_INT    dim_xpw; //[X,P,W]中的向量个数
    GCGE_INT    last_dim_xpw; //上次GCGE迭代中[X,P,W]中的向量个数
    GCGE_INT    last_dim_x; //上次GCGE迭代中X中的向量个数

    GCGE_INT    conv_bs; //全部收敛的批次的特征值个数
    GCGE_INT    unconv_bs; //本批次中未收敛的特征值个数

    //子空间操作用到的工作空间
    void        *dense_matrix; //稠密子空间矩阵

    GCGE_DOUBLE *subspace_matrix; //子空间矩阵
    GCGE_DOUBLE *subspace_evec; //子空间特征向量
    GCGE_DOUBLE *subspace_dtmp; //小规模的工作空间，原来的AA_sub等都放到work_space
    GCGE_INT    *subspace_itmp; 

    GCGE_DOUBLE *eval; //子空间矩阵的特征值

    //正交化用到的工作空间
    GCGE_INT    *orth_ind; //用于正交化过程中存储非零向量的编号
}GCGE_WORKSPACE;

//一些结构体的空间的创建与销毁
void GCGE_WORKSPACE_Create(GCGE_WORKSPACE **workspace);
void GCGE_WORKSPACE_Setup(GCGE_WORKSPACE *workspace, GCGE_PARA *para, GCGE_OPS *ops, void *A);
void GCGE_WORKSPACE_Free(GCGE_WORKSPACE **workspace, GCGE_PARA *para, GCGE_OPS *ops);

#endif
