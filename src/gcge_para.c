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
    (*para)->w_orth_type     = "bgs"; //使用稳定的块正交化方法

     //正交化参数
    (*para)->orth_para = (GCGE_ORTH_PARA *)malloc(sizeof(GCGE_PARA));
    (*para)->orth_para->orth_zero_tol   = 1e-16;
    (*para)->orth_para->reorth_tol      = 0.75;
    (*para)->orth_para->max_reorth_time = 3;
    (*para)->orth_para->criterion_tol   = 1e-6;
    (*para)->orth_para->print_orth_zero = 0;

    (*para)->multi_tol       = 0.2;
    (*para)->if_use_cg       = 1; //使用内置CG迭代计算W
    (*para)->cg_max_it       = 20;
    (*para)->cg_rate         = 1e-2;
     //cg type: 1: 普通的形式, 2: 并行计算的形式(向量内积计算放在一起计算)
    (*para)->cg_type         = 1; 

    (*para)->num_iter        = 0; //output, old:-2
    (*para)->num_unlock      = 6; //output

    (*para)->use_mpi_bcast   = 1; //Rayleigh-Ritz广播子空间特征向量

    (*para)->print_cg_error  = 0;
    (*para)->print_eval      = 1;
    //(*para)->print_matlab    = 0;
    (*para)->print_part_time = 1;
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
    (*para)->dirichlet_boundary  = 0;

    GCGE_STATISTIC_PARA_Create(&((*para)->stat_para));
    return;
}

//把para占用的内存空间释放掉
void GCGE_PARA_Free(GCGE_PARA **para)
{
    GCGE_STATISTIC_PARA_Free(&((*para)->stat_para));
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
        else if(0 == strcmp(argv[arg_index], "-gcge_w_orth_type")) 
        {
            arg_index++;
            para->w_orth_type = argv[arg_index++];
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
       GCGE_Printf("\n");
       GCGE_Printf("Usage: %s [<options>]\n", argv[0]);
       GCGE_Printf("\n");
       GCGE_Printf("  -gcge_nev             <i>  : number of eigenpairs you need              (default: 6)\n");
       GCGE_Printf("  -gcge_ev_max_it       <i>  : maximum of gcg iterations                  (default: 30)\n");
       GCGE_Printf("  -gcge_block_size      <i>  : number of eigenpairs computed in one patch (default: nev)\n");
       GCGE_Printf("  -gcge_if_lobpcg       <i>  : use lobpcg(1) or gcg(0)                    (default: 0)\n");
       GCGE_Printf("  -gcge_given_init_evec <i>  : giben initial eigenvectors or not          (default: 0)\n");
       GCGE_Printf("  -gcge_ev_tol          <d>  : convergence tolerance                      (default: 1e-4)\n");
       GCGE_Printf("  -gcge_conv_type       <c>  : use reletive or abosolute residual         (default: R)\n");
       GCGE_Printf("  -gcge_orth_type       <c>  : use A norm or B norm orthogonal            (default: B)\n");
       GCGE_Printf("  -gcge_w_orth_type     <c>  : use which kind of orthogonalization        (default: bgs)\n");
       GCGE_Printf("  -gcge_orth_zero_tol   <d>  : zero tolerance in orthogonal               (default: 1e-16)\n");
       GCGE_Printf("  -gcge_reorth_tol      <d>  : reorthgonal tolerance                      (default: 0.75)\n");
       GCGE_Printf("  -gcge_max_reorth_time <i>  : maximun reorthogonal times                 (default: 3)\n");
       GCGE_Printf("  -gcge_print_orth_zero <i>  : print the zero index in orthogonal or not  (default: 0)\n");
       GCGE_Printf("  -gcge_multi_tol       <d>  : tolerance for eigenvalue multiplicity      (default: 0.2)\n");
       GCGE_Printf("  -gcge_if_use_cg       <i>  : use the internal cg or not                 (default: 1)\n");
       GCGE_Printf("  -gcge_cg_max_it       <i>  : maximun numbers of cg iterations           (default: 20)\n");
       GCGE_Printf("  -gcge_cg_rate         <d>  : descent rate of residual in cg             (default: 1e-2)\n");
       GCGE_Printf("  -gcge_print_cg_error  <i>  : print residual error in cg or not          (default: 0)\n");
       GCGE_Printf("  -gcge_print_eval      <i>  : print eigenvalue in each iteration or not  (default: 1)\n");
       GCGE_Printf("  -gcge_print_part_time <i>  : print time of each part in gcg or not      (default: 0)\n");
       GCGE_Printf("  -gcge_print_para      <i>  : print the parameters not                   (default: 1)\n");
       GCGE_Printf("  -gcge_print_result    <i>  : print the final result or not              (default: 1)\n");
       GCGE_Printf("\n");
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
//setup感觉应该放在create内，没有做什么特别的事情，并且也不是用户可以直接接触到的程序
void GCGE_PARA_Setup(GCGE_PARA *para)
{
    GCGE_INT nev = para->nev;
     //设置初始没有收敛的特征值的个数，这里设置成nev呢还是再多算一些特征值？ 
     //应该再多算一些特征值    
    para->num_unlock = nev;
     //分批计算特征值的个数: 默认情况下设置为求解特征值的个数
    if(para->block_size == 0)
    {
        para->block_size = nev;
    }
    para->res = (GCGE_DOUBLE*)calloc(nev, sizeof(GCGE_DOUBLE));
    GCGE_INT i = 0;
    //每个特征值相应的残差大小     
    for(i=0; i<nev; i++)
    {
        para->res[i] = 1.0;
    }
}//end for GCGE_PARA_Setup


//打印收敛性信息
void GCGE_PrintIterationInfo(GCGE_DOUBLE *eval, GCGE_PARA *para)
{
    GCGE_INT i = 0;
    GCGE_INT num_iter   = para->num_iter; 
    GCGE_INT num_unlock = para->num_unlock;
    GCGE_STATISTIC_PARA *stat_para = para->stat_para;
    GCGE_Printf("\nnum_unlock(%d)= %d; max_res(%d)= %e; min_res(%d)= %e;\n", 
            num_iter, num_unlock, 
            num_iter, para->max_res,
            num_iter, para->min_res);
    if(para->print_eval == 1)
    {
        GCGE_INT     nev = para->nev; 
        GCGE_DOUBLE  ev_tol = para->ev_tol; 
        GCGE_DOUBLE *res = para->res;
        char        *conv_type;
        if(strcmp(para->conv_type, "R") == 0)
        {
            conv_type = "relative_res";
        }
        else
        {
            conv_type = "abosolute_res";
        }
        for( i=0; i < nev; i++ )
        {
            GCGE_Printf("eval(%d,%3d) = %20.15e ; "
                    "%s(%d,%3d) =  %20.15e ;", 
                    num_iter, i+1, eval[i], 
                    conv_type, num_iter, i+1, res[i]);
            if(res[i] < ev_tol)
            {
                GCGE_Printf(" (converged)\n");
            }
            else if(i == para->max_ind)
            {
                GCGE_Printf(" (max)\n");
            }
            else if(i == para->min_ind)
            {
                GCGE_Printf(" (min)\n");
            }
            else
            {
                GCGE_Printf("\n");
            }
        }
    }
    if(para->print_part_time == 1)
    {
        GCGE_DOUBLE iter_time = 0.0;
        GCGE_PART_TIME *part_time = stat_para->part_time_one_iter;
        iter_time += part_time->x_axpy_time;
        iter_time += part_time->p_axpy_time;
        iter_time += part_time->p_orth_time;
        iter_time += part_time->w_line_time; 
        iter_time += part_time->w_orth_time;
        iter_time += part_time->rr_mat_time;
        iter_time += part_time->rr_eigen_time;
        iter_time += part_time->conv_time;
        GCGE_Printf("\niteration: %d, time: %12.3lf\n",
                num_iter, iter_time);
        GCGE_Printf("\nDetail time and percentage in this iteration:\n"
                "       x_axpy"
                "       p_orth"
                "       p_axpy"
                "    w_lsolver"
                "       w_orth"
                "     rr_eigen"
                "       rr_mat"
                "   check_conv\n"
                " %11.3lf  %11.3lf  %11.3lf  %11.3lf  %11.3lf  %11.3lf  %11.3lf  %11.3lf\n",
                part_time->x_axpy_time,
                part_time->p_orth_time,
                part_time->p_axpy_time,
                part_time->w_line_time, 
                part_time->w_orth_time,
                part_time->rr_eigen_time,
                part_time->rr_mat_time,
                part_time->conv_time);
        if(iter_time < 1e-10)
        {
            GCGE_Printf("Total time in this iteration is less than 1e-10, the percentege is not presented.\n");
        }
        else
        {
            GCGE_Printf(" %11.3lf%% %11.3lf%% %11.3lf%% %11.3lf%% %11.3lf%% %11.3lf%% %11.3lf%% %11.3lf%%\n\n", 
                    100*(part_time->x_axpy_time)  /iter_time,
                    100*(part_time->p_orth_time)  /iter_time,
                    100*(part_time->p_axpy_time)  /iter_time,
                    100*(part_time->w_line_time)  /iter_time,
                    100*(part_time->w_orth_time)  /iter_time,
                    100*(part_time->rr_eigen_time)/iter_time,
                    100*(part_time->rr_mat_time)  /iter_time,
                    100*(part_time->conv_time)    /iter_time);
        }
    }
    fflush(stdout);
}

//GCG迭代结束后打印特征值及收敛情况，时间信息
void GCGE_PrintFinalInfo(GCGE_DOUBLE *eval, GCGE_PARA *para)
{
    GCGE_INT    i = 0;
    GCGE_INT    nev = para->nev;
    GCGE_DOUBLE ev_tol = para->ev_tol;
    GCGE_DOUBLE *res = para->res;
    GCGE_PART_TIME *part_time = para->stat_para->part_time_total;
    GCGE_DOUBLE total_time = 0.0;
    total_time += part_time->x_axpy_time;
    total_time += part_time->p_axpy_time;
    total_time += part_time->p_orth_time;
    total_time += part_time->w_line_time; 
    total_time += part_time->w_orth_time;
    total_time += part_time->rr_mat_time;
    total_time += part_time->rr_eigen_time;
    total_time += part_time->conv_time;

    GCGE_Printf("\nFinal eigenvalues and residual information:\n");
    GCGE_Printf( "\n            eigenvalue           ");
    if(strcmp(para->conv_type, "R") == 0)
    {
        GCGE_Printf("relative residual\n\n");
    }
    else
    {
        GCGE_Printf("abosolute residual\n\n");
    }
    for( i=0; i < nev; i++ )
    {
        GCGE_Printf("%5d  %20.15e  %20.15e ",
                i+1, eval[i], res[i]);
        if(res[i] < ev_tol)
        {
            GCGE_Printf(" (converged)\n");
        }
        else
        {
            GCGE_Printf("\n");
        }
    }
    GCGE_Printf("\nTotal gcg iteration: %d, total time: %.3lf, ",
            para->num_iter, total_time);
    if(para->num_unlock == 0)
    {
        GCGE_Printf("all converged!\n");
    }
    else
    {
        GCGE_Printf("%d converged, %d unconverged.\n",
                para->nev - para->num_unlock,
                para->num_unlock);
    }
    GCGE_Printf("\nDetail time in total gcg iteration:\n"
            "       x_axpy"
            "       p_orth"
            "       p_axpy"
            "    w_lsolver"
            "       w_orth"
            "     rr_eigen"
            "       rr_mat"
            "   check_conv\n"
            " %12.3lf %12.3lf %12.3lf %12.3lf %12.3lf %12.3lf %12.3lf %12.3lf\n",
            part_time->x_axpy_time,
            part_time->p_orth_time,
            part_time->p_axpy_time,
            part_time->w_line_time, 
            part_time->w_orth_time,
            part_time->rr_eigen_time,
            part_time->rr_mat_time,
            part_time->conv_time);
    if(total_time < 1e-10)
    {
        GCGE_Printf("Total time is less than 1e-10, the percentege is not presented.\n");
    }
    else
    {
        GCGE_Printf(" %11.3lf%% %11.3lf%% %11.3lf%% %11.3lf%% %11.3lf%% %11.3lf%% %11.3lf%% %11.3lf%%\n\n", 
                100*(part_time->x_axpy_time)  /total_time,
                100*(part_time->p_orth_time)  /total_time,
                100*(part_time->p_axpy_time)  /total_time,
                100*(part_time->w_line_time)  /total_time,
                100*(part_time->w_orth_time)  /total_time,
                100*(part_time->rr_eigen_time)/total_time,
                100*(part_time->rr_mat_time)  /total_time,
                100*(part_time->conv_time)    /total_time);
    }
}

//打印参数信息
void GCGE_PrintParaInfo(GCGE_PARA *para)
{
    if(para->print_para == 1)
    {
       GCGE_Printf("\n");
       GCGE_Printf("GCGE_PARA:\n");
       GCGE_Printf("\n");
       GCGE_Printf("  nev             : %d, (number of eigenpairs you need              )\n", para->nev            );
       GCGE_Printf("  ev_max_it       : %d, (maximum of gcg iterations                  )\n", para->ev_max_it      );
       GCGE_Printf("  block_size      : %d, (number of eigenpairs computed in one patch )\n", para->block_size     );
       GCGE_Printf("  given_init_evec : %d, (given initial eigenvectors or not          )\n", para->given_init_evec);
       GCGE_Printf("  max_reorth_time : %d, (maximun reorthogonal times                 )\n", para->orth_para->max_reorth_time);
       GCGE_Printf("  print_orth_zero : %d, (print the zero index in orthogonal or not  )\n", para->orth_para->print_orth_zero);
       GCGE_Printf("  if_use_cg       : %d, (use the internal cg or not                 )\n", para->if_use_cg      );
       GCGE_Printf("  cg_max_it       : %d, (maximun numbers of cg iterations           )\n", para->cg_max_it      );
       GCGE_Printf("  print_cg_error  : %d, (print residual error in cg or not          )\n", para->print_cg_error );
       GCGE_Printf("  print_eval      : %d, (print eigenvalue in each iteration or not  )\n", para->print_eval     );
       GCGE_Printf("  print_part_time : %d, (print time of each part in gcg or not      )\n", para->print_part_time);
       GCGE_Printf("  print_para      : %d, (print the parameters not                   )\n", para->print_para     );
       GCGE_Printf("  print_result    : %d, (print the final result or not              )\n", para->print_result   );
       GCGE_Printf("  conv_type       : %s, (use reletive or abosolute residual         )\n", para->conv_type      );
       GCGE_Printf("  orth_type       : %s, (use A norm or B norm orthogonal            )\n", para->orth_type      );
       GCGE_Printf("  w_orth_type     : %s, (use which kind of orthogonalization        )\n", para->w_orth_type    );
       GCGE_Printf("  ev_tol          : %e, (convergence tolerance                      )\n", para->ev_tol         );
       GCGE_Printf("  orth_zero_tol   : %e, (zero tolerance in orthogonal               )\n", para->orth_para->orth_zero_tol  );
       GCGE_Printf("  reorth_tol      : %e, (reorthgonal tolerance                      )\n", para->orth_para->reorth_tol     );
       GCGE_Printf("  cg_rate         : %e, (descent rate of residual in cg             )\n", para->cg_rate        );
       GCGE_Printf("  multi_tol       : %f, (tolerance for eigenvalue multiplicity      )\n", para->multi_tol      );
       GCGE_Printf("\n");
    }
}
