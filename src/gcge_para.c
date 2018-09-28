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
void GCGE_PARA_Create(GCGE_PARA **para)
{
    /* TODO Should class similar para -> ls_para, orth_para, es_para, print_para */
    (*para) = (GCGE_PARA *)malloc(sizeof(GCGE_PARA));
    //默认计算6个特征值
    (*para)->nev             = 1;
    (*para)->ev_max_it       = 30;
    (*para)->block_size      = 0;

    (*para)->if_lobgcg       = 0; //用gcg
    (*para)->given_init_evec = 0;
    (*para)->ev_tol          = 1e-4;

    (*para)->conv_type       = "R"; //使用相对残差判断收敛性
    (*para)->orth_type       = "B"; //使用B正交

    //正交化参数
    (*para)->orth_para = (GCGE_ORTH_PARA *)malloc(sizeof(GCGE_PARA));
    (*para)->orth_para->orth_zero_tol   = 1e-16;
    (*para)->orth_para->reorth_tol      = 0.75;
    (*para)->orth_para->max_reorth_time = 3;
    (*para)->orth_para->print_orth_zero = 0;

    (*para)->multi_tol       = 0.2;
    (*para)->if_use_cg       = 1; //使用内置CG迭代计算W
    (*para)->cg_max_it       = 20;
    (*para)->cg_rate         = 1e-2;

    (*para)->num_iter        = 0; //output, old:-2
    (*para)->num_unlock      = 6; //output

    (*para)->print_cg_error  = 0;
    (*para)->print_eval      = 1;
    //(*para)->print_matlab    = 0;
    (*para)->print_part_time = 0;
    (*para)->print_para      = 1;
    (*para)->print_result    = 1;
    /* print_level: 0：什么都不打印
                    1：打印最终结果
                    2：多打印每次迭代的特征值和残差
                    3：多打印每次迭代的收敛个数，最大最小残差
                    4：多打印每次迭代的各部分时间及占比
                    5：多打印正交化中0向量出现的位置 */
    //(*para)->print_level     = 4; /* TODO */

    (*para)->max_res         = 0.0;
    (*para)->max_ind         = 0;
    (*para)->min_res         = 0.0;
    (*para)->min_ind         = 0;
    (*para)->sum_res         = 0.0;
    return;
}

void GCGE_PARA_Free(GCGE_PARA **para)
{
    free((*para)->orth_para);  (*para)->orth_para = NULL;
    if((*para)->res)
    {
        free((*para)->res);  (*para)->res = NULL;
    }
    free(*para);  *para = NULL;
    return;
}

GCGE_INT GCGE_PARA_SetFromCommandLine(GCGE_PARA *para, GCGE_INT argc, char **argv)
{
    GCGE_INT arg_index = 0, print_usage = 0;
  
    while(arg_index < argc) 
    {
        if(0 == strcmp(argv[arg_index], "-gcge_nev")) 
        {
            arg_index++;
            para->nev = atoi(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_ev_max_it")) 
        {
            arg_index++;
            para->ev_max_it = atoi(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_block_size")) 
        {
            arg_index++;
            para->block_size = atoi(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_if_lobgcg")) 
        {
            arg_index++;
            para->if_lobgcg = atoi(argv[arg_index++]);
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
        else if(0 == strcmp(argv[arg_index], "-gcge_orth_zero_tol")) 
        {
            arg_index++;
            para->orth_para->orth_zero_tol = atof(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_reorth_tol")) 
        {
            arg_index++;
            para->orth_para->reorth_tol = atof(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_max_reorth_time")) 
        {
            arg_index++;
            para->orth_para->max_reorth_time = atoi(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_print_orth_zero")) 
        {
            arg_index++;
            para->orth_para->print_orth_zero = atoi(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_multi_tol")) 
        {
            arg_index++;
            para->multi_tol = atof(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_if_use_cg")) 
        {
            arg_index++;
            para->if_use_cg = atoi(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_cg_max_it")) 
        {
            arg_index++;
            para->cg_max_it = atoi(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_cg_rate")) 
        {
            arg_index++;
            para->cg_rate = atoi(argv[arg_index++]);
        }
        else if(0 == strcmp(argv[arg_index], "-gcge_print_cg_error")) 
        {
            arg_index++;
            para->print_cg_error = atoi(argv[arg_index++]);
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
        else if(0 == strcmp(argv[arg_index], "-gcge_print_para")) 
        {
            arg_index++;
            para->print_para = atoi(argv[arg_index++]);
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
       printf("  -gcge_ev_max_it       <i>  : maximum of gcg iterations                  (default: 30)\n");
       printf("  -gcge_block_size      <i>  : number of eigenpairs computed in one patch (default: nev)\n");
       printf("  -gcge_if_lobpcg       <i>  : use lobpcg(1) or gcg(0)                    (default: 0)\n");
       printf("  -gcge_given_init_evec <i>  : giben initial eigenvectors or not          (default: 0)\n");
       printf("  -gcge_ev_tol          <d>  : convergence tolerance                      (default: 1e-4)\n");
       printf("  -gcge_conv_type       <c>  : use reletive or abosolute residual         (default: R)\n");
       printf("  -gcge_orth_type       <c>  : use A norm or B norm orthogonal            (default: B)\n");
       printf("  -gcge_orth_zero_tol   <d>  : zero tolerance in orthogonal               (default: 1e-16)\n");
       printf("  -gcge_reorth_tol      <d>  : reorthgonal tolerance                      (default: 0.75)\n");
       printf("  -gcge_max_reorth_time <i>  : maximun reorthogonal times                 (default: 3)\n");
       printf("  -gcge_print_orth_zero <i>  : print the zero index in orthogonal or not  (default: 0)\n");
       printf("  -gcge_multi_tol       <d>  : tolerance for eigenvalue multiplicity      (default: 0.2)\n");
       printf("  -gcge_if_use_cg       <i>  : use the internal cg or not                 (default: 1)\n");
       printf("  -gcge_cg_max_it       <i>  : maximun numbers of cg iterations           (default: 20)\n");
       printf("  -gcge_cg_rate         <d>  : descent rate of residual in cg             (default: 1e-2)\n");
       printf("  -gcge_print_cg_error  <i>  : print residual error in cg or not          (default: 0)\n");
       printf("  -gcge_print_eval      <i>  : print eigenvalue in each iteration or not  (default: 1)\n");
       printf("  -gcge_print_part_time <i>  : print time of each part in gcg or not      (default: 0)\n");
       printf("  -gcge_print_para      <i>  : print the parameters not                   (default: 1)\n");
       printf("  -gcge_print_result    <i>  : print the final result or not              (default: 1)\n");
       printf("\n");
    }
    if (print_usage)
       return 1;
    else 
       return 0;
}

void GCGE_PARA_SetNumEigen(GCGE_PARA *para, GCGE_INT nev)
{
    para->nev = nev;
}

void GCGE_PARA_Setup(GCGE_PARA *para)
{
    GCGE_INT nev = para->nev;
    para->num_unlock = nev;
    if(para->block_size == 0)
    {
        para->block_size = nev;
    }
    para->res = (GCGE_DOUBLE*)calloc(nev, sizeof(GCGE_DOUBLE));
    GCGE_INT i = 0;
    for(i=0; i<nev; i++)
    {
        para->res[i] = 1.0;
    }
}

