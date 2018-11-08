/*
 * =====================================================================================
 *
 *       Filename:  test_orthogonal.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/27/2018 04:47:29 PM
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

/* TODO should be modified orthogonal's workspace shoud be VEC** */
#include "gcge_solver.h"

#include "gcge_app_csr.h"

/**
 * @brief 
 *
 * 正交化V, 从start到end, 且返回end, 因为正交化过程中会去掉0向量
 *
 * V_tmp只需要一个，记录当前要正交的向量V[current]与B的乘积
 * d_tmp记录当前要正交的向量V[current]与之前向量以及与自己的B内积, 长度为输入的end
 *
 * 过程中需要做重正交化，当向量的模改变太大(此时，会包含太多以前的分量)
 *
 * 若向量的模太小，则与*end向量交换指针，并(*end)--
 *
 * 最后做一次归一化
 *
 * 注意：大写表示向量组，小写表示单向量
 *
 * for current = start ; current < (*end); ++current
 *     计算V_tmp
 *     reorth_count = 0;
 *     do{
 *     MatDotVec(B, V[current], V_tmp)
 *     计算d_tmp
 *
 *     mv_s = {0, 0}
 *     mv_e = {current+1, 1}
 *     MultiVecInnerProd(V, &V_tmp, d_tmp, "ns", d_tmp, mv_s, mv_e, current+1, ops )
 *
 *    // V[current] = V[current] - sum_{i=0}^{current-1} d_tmp[i] V[i]
 *    // d_tmp[current] = d_tmp[current] - sum_{i=0}^{current-1} d_tmp[i]^2
 *    // V[current] = V[current] / d_tmp[current]
 * 
 *     if(current != 0)
 *     {
 *       mv_s = {0, 0}
 *       mv_e = {current, 1}
 *       MultiVecLinearComb (V, &V_tmp, mv_s, mv_e, d_tmp, current, ops)
 *       VecAxpby(-1, V_tmp, 1, V[current]);
 *     }
 *
 *     vin   = d_tmp[current] 
 *     d_tmp[current] = d_tmp[current] - sum_{i=0}^{current-1} d_tmp[i]^2
 *
 *     ratio = d_tmp[current] / vin
 *
 *     reorth_count++
 *     
 *     }while(ratio < orth_para->reorth_tol && reorth_count < orth_para->max_reorth_count)
 *
 *     if (d_tmp[current] > orth_para->orth_zero_tol)
 *     {
 *         VecAxpby(0, V_tmp, 1/d_tmp[current], V[current])
 *     }
 *     else 
 *     {
 *         mv_s = {current, *end-1}
 *         mv_e = {current+1, *end}
 *         MultiVecSwap(V, V, mv_s, mv_e, ops);
 *         (*end)--;
 *         current--;
 *     }
 * end
 *
 * @param V
 * @param start 从0 到 start-1已经正交归一
 * @param end   输入要正交的向量个数， 返回被正交化的向量个数
 * @param B
 * @param ops
 * @param orth_para
 * @param V_tmp 长度1
 * @param d_tmp 长度×end
 */
void Orthogonal(void **V, GCGE_INT start, GCGE_INT *end, 
      void *B, GCGE_OPS *ops, GCGE_ORTH_PARA *orth_para, void **V_tmp, GCGE_DOUBLE *d_tmp)
{
    GCGE_INT    current = 0;
    GCGE_INT    reorth_count = 0;
    GCGE_INT    mv_s[2] = { 0,0 };
    GCGE_INT    mv_e[2] = { 0,0 };
    GCGE_DOUBLE vin = 0.0;
    GCGE_DOUBLE ratio = 1.0;
    void        *vec_current;
    void        *vec_tmp;
    GCGE_INT    i = 0;
    for(current = start; current < (*end); ++current)
    {
        reorth_count = 0;
        do{
            ops->GetVecFromMultiVec(V_tmp, 0, &vec_tmp);
            ops->GetVecFromMultiVec(V, current, &vec_current);
            if(B == NULL)
            {
                ops->VecAxpby(1.0, vec_current, 0.0, vec_tmp);
            }
            else
            {
                ops->MatDotVec(B, vec_current, vec_tmp);
            }
            ops->RestoreVecForMultiVec(V_tmp, 0, &vec_tmp);
            ops->RestoreVecForMultiVec(V, current, &vec_current);
            
            mv_s[0] = 0;
            mv_s[1] = 0;
            mv_e[0] = current+1;
            mv_e[1] = 1;
            ops->MultiVecInnerProd(V, V_tmp, d_tmp, "ns", mv_s, mv_e, current+1, ops);
            
            if(current != 0)
            {
                mv_s[0] = 0;
                mv_s[1] = 0;
                mv_e[0] = current;
                mv_e[1] = 1;
                ops->MultiVecLinearComb(V, V_tmp, mv_s, mv_e, d_tmp, current, NULL, -1, ops);
                
                ops->GetVecFromMultiVec(V_tmp, 0, &vec_tmp);
                ops->GetVecFromMultiVec(V, current, &vec_current);
                ops->VecAxpby(-1.0, vec_tmp, 1.0, vec_current);
                ops->RestoreVecForMultiVec(V_tmp, 0, &vec_tmp);
                ops->RestoreVecForMultiVec(V, current, &vec_current);
            }
            
            vin = sqrt(d_tmp[current]);
            d_tmp[current] = d_tmp[current] - GCGE_VecDotVecSubspace(d_tmp, d_tmp, current);
            ratio = sqrt(d_tmp[current])/vin;
            //printf("reorth_count: %d, vin, vout: %f, %f\n", reorth_count, vin, sqrt(d_tmp[current]));

            reorth_count++;
            
        }while(ratio < orth_para->reorth_tol && reorth_count < orth_para->max_reorth_time);
 
      if(d_tmp[current] > orth_para->orth_zero_tol)
      {
          ops->GetVecFromMultiVec(V, current, &vec_current);
          ops->VecAxpby(0, vec_current, 1/sqrt(d_tmp[current]), vec_current);
          ops->RestoreVecForMultiVec(V, current, &vec_current);
      }
      else 
      {
          mv_s[0] = current;
          mv_s[1] = *end-1;
          mv_e[0] = current+1;
          mv_e[1] = *end;
          ops->MultiVecSwap(V, V, mv_s, mv_e, ops);
          (*end)--;
          current--;
          if(orth_para->print_orth_zero == 1)
          {
              printf("In Orthogonal, there is a zero vector!, "
                      "current = %d, start = %d, end = %d\n", current, start, *end);
          }
      }
    }
}

int main(int argc, char* argv[])
{
    srand((unsigned)time(NULL));
    GCGE_SOLVER *solver;
    GCGE_INT error = GCGE_SOLVER_Create(&solver, argc, argv);
    if(error)
    {
        GCGE_SOLVER_Free(&solver);
        exit(0);
    }

    //创建矩阵
    const char *file_A = "../data/testA";
    const char *file_B = "../data/testB";
    CSR_MAT *A = CSR_ReadMatFile(file_A);
    CSR_MAT *B = CSR_ReadMatFile(file_B);

    CSR_PrintMat(A);
    CSR_PrintMat(B);

    //setup
    GCGE_SOLVER_SetMatA(solver, A);
    GCGE_SOLVER_SetMatB(solver, B);
    GCGE_SOLVER_SetCSROps(solver);
    GCGE_SOLVER_Setup(solver);

    int num_vec = 3;
    CSR_VEC **multi_vec = (CSR_VEC**)malloc(num_vec*sizeof(CSR_VEC*));
    int i = 0; 
    for(i=0; i<num_vec; i++)
    {
        CSR_BuildVecByMat(A, multi_vec+i);
    }
    CSR_VEC **multi_vec_2 = (CSR_VEC**)malloc(num_vec*sizeof(CSR_VEC*));
    for(i=0; i<num_vec; i++)
    {
        CSR_BuildVecByMat(A, multi_vec_2+i);
    }
    CSR_VEC *vec;
    CSR_BuildVecByMat(A, &vec);

    for(i=0; i<vec->size; i++)
    {
       vec->Entries[i] = i+1;
    }

    solver->ops->MultiVecSetRandomValue(multi_vec, num_vec, solver->ops);
    solver->ops->VecAxpby(1.0, multi_vec[1], 0.0, multi_vec[2]);
    printf("init_vec:\n");
    int k;
    for(k=0; k<num_vec; k++)
    {
        /*
       for(i=0; i<multi_vec[k]->size; i++)
       {
	  multi_vec[k]->Entries[i] = i+k;
       }
       */
       CSR_PrintVec(multi_vec[k]);
    }
    int start = 0, end = 1;
    Orthogonal((void **)multi_vec, start, &end, (void *)A, 
	  solver->ops, solver->para->orth_para, 
          solver->workspace->V_tmp,
          solver->workspace->subspace_dtmp);
    /* SHOULD be modified to ops */
    //GCGE_Orthogonal((void **)multi_vec, start, &end, (void *)B, 
	//  solver->ops, solver->para->orth_para, solver->workspace);
    //Orthogonal((void **)multi_vec, start, &end, NULL, 
    start = 1, end = 3;
    Orthogonal((void **)multi_vec, start, &end, (void *)A, 
	  solver->ops, solver->para->orth_para, 
          solver->workspace->V_tmp,
          solver->workspace->subspace_dtmp);
    printf("end: %d\n", end);

    for(k=0; k<num_vec; k++)
    {
       CSR_PrintVec(multi_vec[k]);
    }

    double vTv;
    solver->ops->VecInnerProd((void *)multi_vec[2], (void *)multi_vec[2], &vTv);
    printf ( "vTv = %f\n", vTv );

    double mvTmv[9];
    int    mv_s[2] = {0, 0};
    int    mv_e[2] = {3, 3};
    solver->ops->MatDotMultiVec((void*)A, (void **)multi_vec, (void **)multi_vec_2, 
            mv_s, mv_e, solver->ops);
    solver->ops->MultiVecInnerProd((void **)multi_vec, (void **)multi_vec_2, (double *)&mvTmv, 
    	  "S", mv_s, mv_e, 3, solver->ops);
    //solver->ops->MultiVecInnerProd((void **)multi_vec, (void **)multi_vec, (double *)&mvTmv, 
    //	  "S", mv_s, mv_e, 3, solver->ops);

    printf ( "mvTmv\n" );
    //printf ( "%f\t%f\n%f\t%f\n", mvTmv[0], mvTmv[1], mvTmv[2], mvTmv[3] );
    printf ( "%f\t%f\t%f\n", mvTmv[0], mvTmv[3], mvTmv[6] );
    printf ( "%f\t%f\t%f\n", mvTmv[1], mvTmv[4], mvTmv[7] );
    printf ( "%f\t%f\t%f\n", mvTmv[2], mvTmv[5], mvTmv[8] );
    
    double entries[6] = {0.0, 1.0, 2.0, 1.0, 2.0, 3.0};
    end = 2;
    GCGE_OrthogonalSubspace((double *)entries, 0, &end, 3, 
	  NULL, -1, solver->para->orth_para, solver->workspace);
    //printf ( "entries\n" );
    //printf ( "%f\t%f\t%f\n%f\t%f\t%f\n", 
    //	  entries[0], entries[1], entries[2], 
    //	  entries[3], entries[4], entries[5]  );


    //释放矩阵空间
    GCGE_SOLVER_Free(&solver);
    CSR_MatFree(&A);
    CSR_MatFree(&B);
    for(i=0; i<num_vec; i++)
    {
        CSR_VecFree(multi_vec+i);
    }
    free(multi_vec); multi_vec = NULL;
    for(i=0; i<num_vec; i++)
    {
        CSR_VecFree(multi_vec_2+i);
    }
    free(multi_vec_2); multi_vec_2 = NULL;
    CSR_VecFree(&vec);

    return 0;
}
