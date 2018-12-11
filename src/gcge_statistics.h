/*
 * =====================================================================================
 *
 *       Filename:  gcge_statistics.h
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
#ifndef _GCGE_STATISTICS_H_
#define _GCGE_STATISTICS_H_

#include "gcge_config.h"

typedef struct GCGE_PART_TIME_ {
    GCGE_DOUBLE x_axpy_time;
    GCGE_DOUBLE x_orth_time; 
    GCGE_DOUBLE p_axpy_time;
    GCGE_DOUBLE p_orth_time;
    GCGE_DOUBLE w_line_time; /* linear solver */
    GCGE_DOUBLE w_orth_time;
    GCGE_DOUBLE rr_mat_time; /* rayleigh_ritz */
    GCGE_DOUBLE rr_eigen_time;
    GCGE_DOUBLE conv_time;
}GCGE_PART_TIME;

typedef struct GCGE_OPS_TIME_COUNT_ {
    GCGE_DOUBLE mv_time;
    GCGE_INT    mv_count;
    GCGE_DOUBLE ip_time; /** inner product */
    GCGE_INT    ip_count;
    GCGE_DOUBLE axpy_time;
    GCGE_INT    axpy_count;
}GCGE_OPS_TIME_COUNT;

typedef struct GCGE_STATISTIC_PARA_ {
	GCGE_PART_TIME      *part_time_total;
	GCGE_PART_TIME      *part_time_one_iter;
	GCGE_OPS_TIME_COUNT *ops_time_count_total;
	GCGE_OPS_TIME_COUNT *ops_time_count_linear_solver;
}GCGE_STATISTIC_PARA;

void GCGE_STATISTIC_PARA_Create(GCGE_STATISTIC_PARA **stat_para);
void GCGE_STATISTIC_PARA_Free  (GCGE_STATISTIC_PARA **stat_para);

void GCGE_PART_TIME_Create(GCGE_PART_TIME **part_time);
void GCGE_PART_TIME_Free  (GCGE_PART_TIME **part_time);

void GCGE_OPS_TIME_COUNT_Create(GCGE_OPS_TIME_COUNT **ops_time_count);
void GCGE_OPS_TIME_COUNT_Free  (GCGE_OPS_TIME_COUNT **ops_time_count);
#endif

