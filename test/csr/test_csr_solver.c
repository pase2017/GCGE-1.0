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
#include "gcge.h"
#include "gcge_app_csr.h"

int main(int argc, char* argv[])
{
    //mwInit();
    srand((unsigned)time(NULL));
   
    //创建矩阵
    //const char *file_A = "../data/testA";
    //const char *file_B = "../data/testB";
    //const char *file_A = "../test/data/Andrews.txt";
    const char *file_A = "../data/A_5.txt";
    const char *file_B = "../data/M_5.txt";
    CSR_MAT *A = CSR_ReadMatFile(file_A);
    CSR_MAT *B = CSR_ReadMatFile(file_B);
    
    int nev = 20;
    //创建一个特征值solver实例
    GCGE_SOLVER *csr_solver = GCGE_CSR_Solver_Init(A, B, nev, argc,  argv);   
     //一些参数的设置
    //csr_solver->para->ev_tol = 1e-8;
    csr_solver->para->dirichlet_boundary = 1;
    csr_solver->para->cg_max_it = 20;
    csr_solver->para->ev_max_it = 150;
    csr_solver->para->cg_type  = 2;
    csr_solver->para->print_part_time = 0;
    csr_solver->para->print_eval = 0;
     
     //求解特征值问题
    GCGE_SOLVER_Solve(csr_solver);  
    //释放矩阵空间
    CSR_MatFree(&A);
    CSR_MatFree(&B);
    GCGE_SOLVER_Free_All(&csr_solver);
    return 0;
}
