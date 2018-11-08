/*
 * =====================================================================================
 *
 *       Filename:  gcge_statistics.c
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
#include "gcge_statistics.h"

//创建统计算法中时间信息的结构体空间
void GCGE_STATISTIC_PARA_Create(GCGE_STATISTIC_PARA **stat_para)
{
    (*stat_para) = (GCGE_STATISTIC_PARA*)malloc(sizeof(GCGE_STATISTIC_PARA));
    GCGE_PART_TIME_Create(&((*stat_para)->part_time_total));
    GCGE_PART_TIME_Create(&((*stat_para)->part_time_one_iter));
    GCGE_OPS_TIME_COUNT_Create(&((*stat_para)->ops_time_count_total));
    GCGE_OPS_TIME_COUNT_Create(&((*stat_para)->ops_time_count_linear_solver));
}

//释放统计算法中时间信息的结构体空间
void GCGE_STATISTIC_PARA_Free(GCGE_STATISTIC_PARA **stat_para)
{
    GCGE_PART_TIME_Free(&((*stat_para)->part_time_total));
    GCGE_PART_TIME_Free(&((*stat_para)->part_time_one_iter));
    GCGE_OPS_TIME_COUNT_Free(&((*stat_para)->ops_time_count_total));
    GCGE_OPS_TIME_COUNT_Free(&((*stat_para)->ops_time_count_linear_solver));
    free(*stat_para); *stat_para = NULL;
}

//统计GCGE算法各部分时间信息的初始化
void GCGE_PART_TIME_Create(GCGE_PART_TIME **part_time)
{
    (*part_time) = (GCGE_PART_TIME*)malloc(sizeof(GCGE_PART_TIME));
    (*part_time)->w_line_time = 0.0;
    (*part_time)->x_axpy_time = 0.0;
    (*part_time)->x_orth_time = 0.0;
    (*part_time)->p_orth_time = 0.0;
    (*part_time)->p_axpy_time = 0.0;
    (*part_time)->w_orth_time = 0.0;
    (*part_time)->rr_eigen_time = 0.0;
    (*part_time)->rr_mat_time = 0.0;
    (*part_time)->conv_time = 0.0;
}

void GCGE_PART_TIME_Free(GCGE_PART_TIME **part_time)
{
    free(*part_time);  *part_time = NULL;
}

//统计GCGE各种操作时间信息的初始化
void GCGE_OPS_TIME_COUNT_Create(GCGE_OPS_TIME_COUNT **ops_time_count)
{
    (*ops_time_count) = (GCGE_OPS_TIME_COUNT*)malloc(sizeof(GCGE_OPS_TIME_COUNT));
    (*ops_time_count)->mv_time  = 0.0;
    (*ops_time_count)->mv_count = 0;
    (*ops_time_count)->ip_time  = 0.0;
    (*ops_time_count)->ip_count = 0;
    (*ops_time_count)->axpy_time  = 0.0;
    (*ops_time_count)->axpy_count = 0;
}

void GCGE_OPS_TIME_COUNT_Free(GCGE_OPS_TIME_COUNT **ops_time_count)
{
    free(*ops_time_count); *ops_time_count = NULL;
}
