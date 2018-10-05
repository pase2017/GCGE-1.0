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

//GCGE_INPUT_PARA包含用户可以调节的算法相关的参数
//ev表示特征对，eval表示特征值，evec表示特征向量
typedef struct GCGE_EIGSOL_PARA_ {

    GCGE_INT    nev;        //要求解的特征值个数，默认值为６
    GCGE_INT    ev_max_iter;//gcg最大迭代次数
    GCGE_INT    block_size; //分批计算，每批计算的特征值个数, 默认为nev
    GCGE_INT    given_init_evec; //用户是否给定初始特征向量,0:用户不给定,1:用户给定

    GCGE_DOUBLE ev_tol;    //判定特征值已收敛的阈值
    char*       conv_type; //判定收敛的方式,1:绝对残差, 2:相对残差
    char*       orth_type; //确定是用Ａ正交还是Ｂ正交,1:A_ORTH, 2:B_ORTH
                           //orth_type不能放到ORTH_PARA当中，因为是eigsol中要用到的

    GCGE_DOUBLE eval_multi_tol; //default:0.2,判断特征值重数（或距离较近）的阈值
    GCGE_DOUBLE eval_multi_max; //检查特征值重数时最大检测到nev的多少倍

    GCGE_INT    print_eval;   //是否每次GCGE迭代后打印特征值
    GCGE_INT    print_result; //是否打印最终结果以及每次迭代的收敛个数等
    GCGE_INT    print_part_time; //是否打印各部分的时间

}GCGE_EIGSOL_PARA;

GCGE_INT GCGE_EIGSOL_PARA_SetFromCommandLine(GCGE_PARA *para, GCGE_INT argc, char **argv);
//para不需要创建和销毁，只需要初始化
void GCGE_EIGSOL_PARA_Initialize(GCGE_EIGSOL_PARA *para); //初始化参数
void GCGE_EIGSOL_PARA_SetNumEigen(GCGE_EIGSOL_PARA *para, GCGE_INT nev);
void GCGE_EIGSOL_PARA_SetBlockSize(GCGE_EIGSOL_PARA *para, GCGE_INT block_size);
void GCGE_EIGSOL_PARA_SetEigenvalueTolerance(GCGE_EIGSOL_PARA *para, GCGE_DOUBLE ev_tol);
void GCGE_EIGSOL_PARA_SetConvergenceType(GCGE_EIGSOL_PARA *para, GCGE_CHAR *conv_type);
void GCGE_EIGSOL_PARA_SetOrthogonalizeType(GCGE_EIGSOL_PARA *para, GCGE_CHAR *orth_type);
void GCGE_EIGSOL_PARA_SetEvalMultiplicityTolerance(GCGE_EIGSOL_PARA *para, GCGE_DOUBLE *eval_multi_tol);
void GCGE_EIGSOL_PARA_SetEvalMultiplicityMaximumIndex(GCGE_EIGSOL_PARA *para, GCGE_DOUBLE *eval_multi_max);
void GCGE_EIGSOL_PARA_SetPrintEval(GCGE_EIGSOL_PARA *para, GCGE_INT *print_eval);
void GCGE_EIGSOL_PARA_SetPrintResult(GCGE_EIGSOL_PARA *para, GCGE_INT *print_result);
void GCGE_EIGSOL_PARA_SetPrintPartTime(GCGE_EIGSOL_PARA *para, GCGE_INT *print_part_time);

//为了防止参数互斥，比如我们设置nev的时候，就把可能会互斥的参数也设置一遍
//对于nev会影响的参数block_size，设置的时候，如果互斥了，要给一个提示，block_size不能比nev大
//要考虑到参数设置的先后关系，互相影响关系
#endif

