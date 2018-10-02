/*
 * =====================================================================================
 *
 *       Filename:  gcge_para.h
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
 /*  */
#ifndef  _GCGE_PARA_H_
#define  _GCGE_PARA_H_

#include "gcge_type.h"

typedef struct GCGE_ORTH_PARA_ {
    //正交化用到的参数
    GCGE_DOUBLE orth_zero_tol; //确定正交化中出现０向量的阈值
    GCGE_DOUBLE reorth_tol; //确定是否进行重正交化的阈值
    GCGE_INT    max_reorth_time; //最大的重正交化次数
    GCGE_INT    print_orth_zero; //print_orthzero确定是否打印正交化过程中出现0向量的信息,0:不打印，1:打印

}GCGE_ORTH_PARA;

typedef struct GCGE_PARA_ {
    GCGE_INT    nev; //要求解的特征值个数，默认值为６
    GCGE_INT    ev_max_it; //gcg最大迭代次数
    GCGE_INT    block_size; //分批计算，每批计算的特征值个数

    //基本参数
    GCGE_INT    if_lobgcg; //使用gcg或lobgcg,0:gcg,1:lobgcg
    GCGE_INT    given_init_evec; //用户是否给定初始特征向量,0:用户不给定,1:用户给定

    //检查收敛性用到的参数
    GCGE_DOUBLE ev_tol; //判定特征值已收敛的阈值
    char*       conv_type; //判定收敛的方式,1:绝对残差, 2:相对残差
    GCGE_INT    num_unlock; //nunlock表示未锁定的特征对个数(0<=nunlock<=nev)
    //GCGE_INT    unconv_bs; //本次迭代中需要计算 P W 的个数

    //正交化用到的参数
    char*       orth_type; //确定是用Ａ正交还是Ｂ正交,1:A_ORTH, 2:B_ORTH
    GCGE_ORTH_PARA *orth_para;

    GCGE_DOUBLE multi_tol;//default:0.2,判断特征值重数（或距离较近）的阈值

    //GCGE内置LinearSolver(CG)参数
    GCGE_INT    if_use_cg; //是否使用内置CG
    GCGE_INT    cg_max_it; //CG最大迭代次数
    GCGE_DOUBLE cg_rate; //终止CG迭代残差下降比例

    //GCGE迭代次数,因为统计的是[X,P,W]都参与计算的迭代次数，所以从-1开始
    GCGE_INT    num_iter;//初始化为-2
    GCGE_DOUBLE max_res; //最大残差
    GCGE_DOUBLE max_ind; //最大残差的位置
    GCGE_DOUBLE min_res; //最小残差
    GCGE_DOUBLE min_ind; //最小残差的位置
    GCGE_DOUBLE sum_res; //残差和
    GCGE_DOUBLE *res; //残差


    //判定是否打印的参数
    GCGE_INT    print_eval; //是否每次GCGE迭代后打印特征值
    GCGE_INT    print_cg_error; //是否打印每次CG迭代的残差
    GCGE_INT    print_orthzero; //print_orthzero确定是否打印正交化过程中出现0向量的信息,0:不打印，1:打印
    GCGE_INT    print_para; //是否打印GCGE_PARA的参数信息
    GCGE_INT    print_result; //是否打印最终结果以及每次迭代的收敛个数等
    GCGE_INT    print_matlab; //如果要用matlab画图，那么就不要打印converged这些
    GCGE_INT    print_part_time; //是否打印各部分的时间
    GCGE_INT    print_level; //用打印层级来控制打印多少东西
    //GCGE_STATISTIC_PARA *stat_para;

}GCGE_PARA;

/**
 * @brief 
 *
 * @param para
 * @param nev
 */
GCGE_INT GCGE_PARA_SetFromCommandLine(GCGE_PARA *para, GCGE_INT argc, char **argv);
void GCGE_PARA_Create(GCGE_PARA **para);
void GCGE_PARA_Free(GCGE_PARA **para);
void GCGE_PARA_SetNumEigen(GCGE_PARA *para, GCGE_INT nev);
void GCGE_PARA_Setup(GCGE_PARA *para);
#endif

