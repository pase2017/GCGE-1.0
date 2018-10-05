/*
 * =====================================================================================
 *
 *       Filename:  gcge_orthogonal.h
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
#ifndef _GCGE_ORTHOGONAL_H_
#define _GCGE_ORTHOGONAL_H_


#include "gcge_type.h"

#include "gcge_para.h"
#include "gcge_ops.h"

#include "gcge_workspace.h"

typedef struct GCGE_ORTH_PARA_ {

    GCGE_DOUBLE orth_zero_tol; //确定正交化中出现０向量的阈值
    GCGE_DOUBLE reorth_tol; //确定是否进行重正交化的阈值
    GCGE_DOUBLE criterion_tol; //确定是否进行块正交化，当sum_res小于这个值的时候进行单个向量正交化
    GCGE_INT    max_reorth_time; //最大的重正交化次数
    GCGE_INT    print_orth_zero; //print_orthzero确定是否打印正交化过程中出现0向量的信息,0:不打印，1:打印

}GCGE_ORTH_PARA;

//GCGE算法中用到的正交化操作
void GCGE_ClassicalOrthogonalize(void **V, GCGE_INT start, GCGE_INT *end, 
                                 void *B, GCGE_OPS *ops);
void GCGE_ModifiedOrthogonalize(void **V, GCGE_INT start, GCGE_INT *end, 
                                void *B, GCGE_OPS *ops);
void GCGE_BlockedOrthogonalize(void **V, GCGE_INT start, GCGE_INT *end, 
                               void *B, GCGE_OPS *ops);

//子空间正交化也要用到orth_para，所以把orth_para的定义放在这里
void GCGE_ModifiedOrthogonalInSubspace(double *V, GCGE_INT ldV, GCGE_INT start, GCGE_INT *end, 
                                       void *B, GCGE_INT ldB, GCGE_ORTH_PARA *orth_para);
#endif
