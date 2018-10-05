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

typedef struct GCGE_EIGSOL_WS_ {

    void        **V; //V用于存储[X,P,W]基底
    void        **V_tmp; //V_tmp用作临时的向量存储空间
    void        **RitzVec; //RitzVec用于存储Ritz向量
    GCGE_INT    space_size_V; //给V分配空间的大小，向量个数
    GCGE_INT    space_size_V_tmp; //给V_tmp分配空间的大小，向量个数
    GCGE_INT    space_size_RitzVec; //给RitzVec分配空间的大小，向量个数

    GCGE_INT    *unlock; //unlock用于存储未锁定的特征对编号(0-nev)

    GCGE_DOUBLE *subspace_matrix; //子空间矩阵
    GCGE_DOUBLE *subspace_evec; //子空间特征向量
    GCGE_DOUBLE *subspace_eval; //子空间矩阵的特征值
    GCGE_DOUBLE *subspace_dtmp; //小规模的double型工作空间
    GCGE_INT    *subspace_itmp; //小规模的int型工作空间

}GCGE_EIGSOL_WS;
//eigsol工作空间的创建与销毁
void GCGE_EIGSOL_WS_Create(GCGE_EIGSOL_WS **eigsol_ws, GCGE_PARA *para, GCGE_OPS *ops, void *A);
void GCGE_EIGSOL_WS_Free(GCGE_EIGSOL_WS **eigsol_ws, GCGE_OPS *ops);

#endif
