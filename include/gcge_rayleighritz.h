/*
 * =====================================================================================
 *
 *       Filename:  gcge_rayleighritz.h
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
#ifndef _GCGE_RAYLEIGHRITZ_H_
#define _GCGE_RAYLEIGHRITZ_H_

#include <string.h>
#include "gcge_type.h"

#include "gcge_para.h"
#include "gcge_ops.h"

#include "gcge_workspace.h"
#include "gcge_utilities.h"
#include "gcge_matvec.h"
#include "gcge_orthogonal.h"
#if GCGE_USE_MPI
#include <mpi.h>
#endif

//计算 work_space = V^T * A * W
void GCGE_ComputeSubspaceMatrixVTAW(void **V, void *A, GCGE_INT start_W, GCGE_INT end_V, GCGE_DOUBLE *subspace_mat, GCGE_OPS *ops, void **workspace);
//计算 subspace_matrix = V^T * A * V
void GCGE_ComputeSubspaceMatrix(void *A, void **V, 
        GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace);
//计算子空间矩阵特征值 AA * u = \lambda * u
void GCGE_ComputeSubspaceEigenpairs(GCGE_DOUBLE *subspace_matrix, 
        GCGE_DOUBLE *eval, GCGE_DOUBLE *subspace_evec, 
        GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace);

//对特征值和特征向量进行排序
void GCGE_SortEigenpairs(GCGE_DOUBLE *eval, GCGE_DOUBLE *evec, GCGE_INT nev, GCGE_INT ldv, GCGE_DOUBLE *work);
GCGE_DOUBLE GCGE_ModuleMaxDouble(GCGE_DOUBLE *a, GCGE_INT n);
#endif
