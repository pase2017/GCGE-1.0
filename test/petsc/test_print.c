/*
 * =====================================================================================
 *
 *       Filename:  test_solver.c
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
#include <math.h>
#include <time.h>

//#include "memwatch.h"
#include <petscsys.h>
#include <petscviewer.h>
#include <petscmat.h>

#include "gcge.h"
#include "gcge_app_petsc.h"
#include "external_petsc.h"
#include <mpi.h>

static char help[] = "Use GCGE-PETSc to solve an eigensystem Ax=kBx with the matrixes loaded from files.\n";

#if 0
void GCGE_Printf(char *fmt, ...)
{
#if GCGE_USE_MPI
    GCGE_INT myrank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if(0 == myrank) {
        va_list vp;
        va_start(vp, fmt);
        vprintf(fmt, vp);
        va_end(vp);
    }
#else
    va_list vp;
    va_start(vp, fmt);
    vprintf(fmt, vp);
    va_end(vp);
#endif
}


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
#endif

int main(int argc, char* argv[])
{
    //创建矩阵
    Mat A, B;
    PetscErrorCode ierr;

    PetscInitialize(&argc,&argv,(char*)0,help);
    GCGE_PARA *para;
    GCGE_PARA_Create(&para);

    int nev = 30;
    para->nev = 30;
    GCGE_PARA_Setup(para);

    GCGE_DOUBLE *eval = (GCGE_DOUBLE*)calloc(nev, sizeof(GCGE_DOUBLE));
    GCGE_PrintIterationInfo(eval, para);
    GCGE_PrintFinalInfo(eval, para);

    //分别对 GCGE_USE_MPI = 0,1 作了测试
    //printf("printf, test_petsc_solver.c nev: %d\n", nev);
    //GCGE_Printf("GCGE_Printf, test_petsc_solver.c nev: %d\n", nev);

    free(eval); eval = NULL;
    GCGE_PARA_Free(&para);
    ierr = PetscFinalize();

    return ierr;
}
