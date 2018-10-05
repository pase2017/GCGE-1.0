/*
 * =====================================================================================
 *
 *       Filename:  gcge_para.c
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
#include "gcge_para.h"



//创建GCGE_Para结构并进行初始化
void GCGE_EIGSOL_PARA_Initialize(GCGE_EIGSOL_PARA *para)
{
    //默认计算6个特征值
    para->nev             = 1;
    para->ev_max_iter     = 30;
    para->block_size      = 0;

    para->given_init_evec = 0;
    para->ev_tol          = 1e-4;

    para->conv_type       = "R"; //使用相对残差判断收敛性
    para->orth_type       = "B"; //使用B正交

    para->eval_multi_tol  = 0.2;
    para->eval_multi_max  = 0.2;

    para->print_eval      = 1;
    para->print_part_time = 0;
    para->print_result    = 1;

    return;
}

GCGE_INT GCGE_EIGSOL_PARA_SetFromCommandLine(GCGE_EIGSOL_PARA *para, GCGE_INT argc, char **argv)
{
    GCGE_INT arg_index = 0, print_usage = 0;
  
    while(arg_index < argc) 
    {
        if(0 == strcmp(argv[arg_index], "-gcge_nev")) 
        {
            arg_index++;
            para->nev = atoi(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_ev_max_iter")) 
        {
            arg_index++;
            para->ev_max_iter = atoi(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_block_size")) 
        {
            arg_index++;
            para->block_size = atoi(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_given_init_evec")) 
        {
            arg_index++;
            para->given_init_evec = atoi(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_ev_tol")) 
        {
            arg_index++;
            para->ev_tol = atof(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_conv_type")) 
        {
            arg_index++;
            para->conv_type = argv[arg_index++];
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_orth_type")) 
        {
            arg_index++;
            para->orth_type = argv[arg_index++];
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_eval_multi_tol")) 
        {
            arg_index++;
            para->eval_multi_tol = atof(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_eval_multi_max")) 
        {
            arg_index++;
            para->eval_multi_max = atof(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_print_eval")) 
        {
            arg_index++;
            para->print_eval = atoi(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_print_part_time")) 
        {
            arg_index++;
            para->print_part_time = atoi(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_print_result")) 
        {
            arg_index++;
            para->print_result = atoi(argv[arg_index++]);
        }
        else if (0 == strcmp(argv[arg_index], "-help"))
        {
           print_usage = 1;
           break;
        }
        else
        {
            arg_index++;
        }
    }
    if (print_usage)
    {
       printf("\n");
       printf("Usage: %s [<options>]\n", argv[0]);
       printf("\n");
       printf("  -gcge_nev             <i>  : number of eigenpairs you need              (default: 6)\n");
       printf("  -gcge_ev_max_iter     <i>  : maximum of gcg iterations                  (default: 30)\n");
       printf("  -gcge_block_size      <i>  : number of eigenpairs computed in one patch (default: nev)\n");
       printf("  -gcge_given_init_evec <i>  : giben initial eigenvectors or not          (default: 0)\n");
       printf("  -gcge_ev_tol          <d>  : convergence tolerance                      (default: 1e-4)\n");
       printf("  -gcge_conv_type       <c>  : use reletive or abosolute residual         (default: R)\n");
       printf("  -gcge_orth_type       <c>  : use A norm or B norm orthogonal            (default: B)\n");
       printf("  -gcge_eval_multi_tol  <d>  : tolerance for eigenvalue multiplicity      (default: 0.2)\n");
       printf("  -gcge_eval_multi_max  <d>  : maximum rate for eigenvalue multiplicity   (default: 0.2)\n");
       printf("  -gcge_print_eval      <i>  : print eigenvalue in each iteration or not  (default: 1)\n");
       printf("  -gcge_print_part_time <i>  : print time of each part in gcg or not      (default: 0)\n");
       printf("  -gcge_print_result    <i>  : print the final result or not              (default: 1)\n");
       printf("\n");
    }
    if (print_usage)
       return 1;
    else 
       return 0;
}

void GCGE_EIGSOL_PARA_SetNumEigen(GCGE_EIGSOL_PARA *para, GCGE_INT nev)
{
    para->nev = nev;
    //设置特征值个数，同时检查block_size不能大于nev
    if(para->block_size > nev)
    {
        para->block_size = nev;
    }
}

void GCGE_EIGSOL_PARA_SetEvMaxIter(GCGE_EIGSOL_PARA *para, GCGE_INT ev_max_iter)
{
    para->ev_max_iter = ev_max_iter;
}
/*
void GCGE_EIGSOL_PARA_SetBlockSize(GCGE_EIGSOL_PARA *para, GCGE_INT block_size)
{
    //设置block_size，不能大于nev
    if(block_size <= para->nev)
    {
        para->block_size = block_size;
    }
    else
    {
        para->block_size = para->nev;
        printf("The block_size cannot be bigger than nev, it has been set as %d.",
                para->nev);
    }
}
*/
void GCGE_EIGSOL_PARA_SetBlockSize(GCGE_EIGSOL_PARA *para, GCGE_INT block_size)
{
    //不做判断，在使用时再做判断
    para->block_size = block_size;
}

void GCGE_EIGSOL_PARA_SetEigenvalueTolerance(GCGE_EIGSOL_PARA *para, GCGE_DOUBLE ev_tol)
{
    para->ev_tol = ev_tol;
}

void GCGE_EIGSOL_PARA_SetConvergenceType(GCGE_EIGSOL_PARA *para, char *conv_type)
{
    //conv_type,相对残差为"R",其余为绝对残差
    para->conv_type = conv_type;
}

void GCGE_EIGSOL_PARA_SetOrthogonalizeType(GCGE_EIGSOL_PARA *para, char *orth_type)
{
    //正交化类型为"A"或"B"
    para->orth_type = orth_type;
}
void GCGE_EIGSOL_PARA_SetEvalMultiplicityTolerance(GCGE_EIGSOL_PARA *para, GCGE_DOUBLE eval_multi_tol)
{
    para->eval_multi_tol = eval_multi_tol;
}
void GCGE_EIGSOL_PARA_SetEvalMultiplicityMaximumIndex(GCGE_EIGSOL_PARA *para, GCGE_DOUBLE eval_multi_max)
{
    para->eval_multi_max = eval_multi_max;
}
void GCGE_EIGSOL_PARA_SetPrintEval(GCGE_EIGSOL_PARA *para, GCGE_INT print_eval)
{
    para->print_eval;
}
void GCGE_EIGSOL_PARA_SetPrintResult(GCGE_EIGSOL_PARA *para, GCGE_INT print_result)
{
    para->print_result = print_result;
}
void GCGE_EIGSOL_PARA_SetPrintPartTime(GCGE_EIGSOL_PARA *para, GCGE_INT print_part_time)
{
    para->print_part_time = print_part_time;
}
