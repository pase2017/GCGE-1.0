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
 *
 * for current = start ; current < (*end); ++current
 *     计算V_tmp
 *     reorth_count = 0;
 *     do{
 *
 *    // V_tmp = B * V[current]
 *    // d_tmp = V^T * V_tmp
 *    // V[current] = V[current] - sum_{i=0}^{current-1} d_tmp[i] V[i]
 *    // d_tmp[current] = d_tmp[current] - sum_{i=0}^{current-1} d_tmp[i]^2
 *    // V[current] = V[current] / d_tmp[current]
 *
 *     ratio = d_tmp[current] / vin
 *
 *     reorth_count++
 *     
       }while(ratio < orth_para->reorth_tol && reorth_count < orth_para->max_reorth_time && vout > orth_para->orth_zero_tol);
 *
 *     if (d_tmp[current] > orth_para->orth_zero_tol)
 *         V[current] = V[current]/vout
 *     else 
 *         V[current] = V[*end-1]
 *         (*end)--;
 *         current--;
 *     end
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


 /**
  * @brief 
  *
  * 正交化V, 从start到end, 且返回end, 因为正交化过程中会去掉0向量
  * 不妨将 V 记为 [V1, V2], 其中 V1 部分已经是正交归一的
  *
  * 由于计算中可能由于机器误差的累积造成正交化缺失，
  * 为了提高稳定性，整个部分做两次（相当于做了一次重正交化）
  *
  * 需要的计算空间 V_tmp = max(end-start, start) 个向量空间
  *                d_tmp = V1^T * V_tmp (start*(end-start)长度)
  * 将正交化分两部分完成
  *
  * 第一部分：去掉V2中V1方向的分量：V2 = V2 - V1 * ( V1^T * B * V2 )
  *
  *   > V_tmp = B * V2 : 
  *     mv_s = { start, 0 }
  *     mv_e = { end, end-start }
  *     ops->MatDotMultiVec(B, V, V_tmp, mv_s, mv_e, ops)
  *
  *     ^ 如果B == NULL
  *       V_tmp = V2 : 
  *       mv_s = { start, 0 }
  *       mv_e = { end, end-start }
  *       ops->MatDotMultiVec(B, V, V_tmp, mv_s, mv_e, ops)
  *
  *   > d_tmp = V1^T * V_tmp:
  *     mv_s = { 0, 0 }
  *     mv_e = { start, end-start }
  *     lda = start (V1的向量个数)
  *     ops->MultiVecInnerProd(V, V, d_tmp, "nonsym", mv_s, mv_e, lda, ops);
  *
  *   > V_tmp = V1 * d_tmp
  *     mv_s = { 0, 0 }
  *     mv_e = { start, end-start }
  *     lda = start (V1的向量个数)
  *     ops->MultiVecLinearComb(V, V_tmp, mv_s, mv_e, d_tmp, lda, ops)
  *
  *   > V2 = V2 - V_tmp
  *     mv_s = { 0, start }
  *     mv_e = { end-start, end }
  *     ops->MultiVecAxpby(-1.0, V_tmp, 1.0, V, mv_s, mv_e, ops)
  *   
  * 第二部分：对V2自身进行正交化
  *     计算 M = V2^T * B * V2
  *     求解特征值问题 M y = theta y
  *     计算 y_i = 1/sqrt(y_i) * y_i
  *     计算 V2 = V2 * Y
  *
  *   > V_tmp = B * V2 : 
  *     mv_s = { start, 0 }
  *     mv_e = { end, end-start }
  *     ops->MatDotMultiVec(B, V, V_tmp, mv_s, mv_e, ops)
  *
  *     ^ 如果B == NULL
  *       V_tmp = V2 : 
  *       mv_s = { start, 0 }
  *       mv_e = { end, end-start }
  *       ops->MatDotMultiVec(B, V, V_tmp, mv_s, mv_e, ops)
  *
  *     d_tmp = V2^T * V_tmp:
  *     mv_s = { start, 0 }
  *     mv_e = { end, end-start }
  *     lda = end-start (V1的向量个数)
  *     ops->MultiVecInnerProd(V, V_tmp, d_tmp, "sym", mv_s, mv_e, lda, ops);
  *
  *     y,theta: d_mat的输出即为特征向量
  *     dsyev(d_mat, theta, y)
  *
  *     计算 y[i] = 1/sqrt(theta[i]) * y[i]
  *     lda = start 
  *     for(current=start; current<(*end); ++current)
  *       norm = sqrt(theta[current])
  *       if theta[current] > orth_zero_tol
  *           GCGE_ArrayScaleInSubspace(1.0/norm, d_mat[current], lda);
  *       else
  *           GCGE_ArrayCopyInSubspace(d_mat[end], d_mat[current], lda);
  *           current--;
  *           (*end)--;
  *       end
  *     end
  *
  *     计算 V_tmp = V2 * Y
  *     mv_s = { start, 0}
  *     mv_e = { end, end-start }
  *     lda = start-end (V2的向量个数)
  *     ops->MultiVecLinearComb(V, V_tmp, mv_s, mv_e, d_tmp, lda, ops)
  *
  *     V2 = V_tmp
  *     mv_s = { 0, start }
  *     mv_e = { end-start, end }
  *     ops->MultiVecAxpby(1.0, V_tmp, 0.0, V, mv_s, mv_e, ops)
  *
  * @param V
  * @param ldV
  * @param start
  * @param end
  * @param B
  * @param ldB
  * @param ops
  * @param orth_para
  * @param V_tmp
  * @param d_tmp
  */
void OrthonormalizationSubspace(double *V, GCGE_INT ldV, GCGE_INT start, GCGE_INT *end, 
      void *B, GCGE_INT ldB, GCGE_ORTH_PARA *orth_para)
{
    GCGE_INT    current = 0;
    GCGE_INT    reorth_count = 0;
    GCGE_DOUBLE vin = 0.0;
    GCGE_DOUBLE vout = 0.0;
    GCGE_DOUBLE ratio = 1.0;
    GCGE_INT    idx = 0;
    GCGE_DOUBLE ip = 1.0;
    for(current = start; current < (*end); ++current)
    {
        vout = GCGE_ArrayNormInSubspace(V+current*ldV, ldV);
        if(current != 0)
        {
            reorth_count = 0;
            do{
                vin = vout;
                for(idx = 0; idx < current; idx++)
                {
                    ip = GCGE_VecDotVecSubspace(V+idx*ldV, V+current*ldV, ldV);
                    GCGE_VecAXPBYSubspace(-ip, V+idx*ldV, 1.0, V+current*ldV, ldV);
                }
                vout = GCGE_ArrayNormInSubspace(V+current*ldV, ldV);

                ratio = vout/vin;

                reorth_count++;

            }while(ratio < orth_para->reorth_tol && reorth_count < orth_para->max_reorth_time && vout > orth_para->orth_zero_tol);
        }

        if(vout > orth_para->orth_zero_tol)
        {
            GCGE_ArrayScaleInSubspace(1.0/vout, V+current*ldV, ldV);
        }
        else 
        {
            GCGE_ArrayCopyInSubspace(V+(*end -1)*ldV, V+current*ldV, ldV);
            (*end)--;
            current--;
            if(orth_para->print_orth_zero == 1)
            {
                printf("In OrthonormalizationSubspace, there is a zero vector!, "
                        "current = %d, start = %d, end = %d\n", current, start, *end);
            }
        }
    }
}


int main(int argc, char* argv[])
{
    printf("111111\n");
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


    //setup
    GCGE_SOLVER_SetMatA(solver, A);
    GCGE_SOLVER_SetMatB(solver, B);
    GCGE_SOLVER_SetCSROps(solver);
    GCGE_SOLVER_Setup(solver);

    int n = 3;
    double entries[9];
    double ip[9];
    int ldV = n;
    int end = n;
    int start = 0;
    int j = 0;
    int i = 0;
    for( j=0; j<n*n; j++ )
    {
       entries[j] = (double)rand()/((double)RAND_MAX+1);
    }
    printf ( "entries\n" );
    printf ( "%f\t%f\t%f\n"
             "%f\t%f\t%f\n" 
             "%f\t%f\t%f\n", 
    	  entries[0], entries[1], entries[2], 
    	  entries[3], entries[4], entries[5],
    	  entries[6], entries[7], entries[8]  );
    for(i=0; i<n; i++)
    {
        //entries[3+i] = entries[6+i];
        entries[3+i] = entries[i];
    }

    //GCGE_OrthonormalizationSubspace((double *)entries, 0, &end, 3, 
    //	  NULL, -1, solver->para->orth_para, solver->workspace);
    end = 1;
    start = 0;
    OrthonormalizationSubspace(entries, ldV, start, &end, 
      NULL, -1, solver->para->orth_para);
    end = 3;
    start = 1;
    OrthonormalizationSubspace(entries, ldV, start, &end, 
      NULL, -1, solver->para->orth_para);
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            ip[i*n+j] = GCGE_VecDotVecSubspace(entries+i*n, entries+j*n, n);
        }
    }

    printf ( "entries\n" );
    printf ( "%f\t%f\t%f\n"
             "%f\t%f\t%f\n" 
             "%f\t%f\t%f\n", 
    	  entries[0], entries[1], entries[2], 
    	  entries[3], entries[4], entries[5],
    	  entries[6], entries[7], entries[8]  );
    printf ( "ip\n" );
    printf ( "%f\t%f\t%f\n"
             "%f\t%f\t%f\n" 
             "%f\t%f\t%f\n", 
    	  ip[0], ip[1], ip[2], 
    	  ip[3], ip[4], ip[5],
    	  ip[6], ip[7], ip[8]  );

    //释放矩阵空间
    GCGE_SOLVER_Free(&solver);
    return 0;
}
