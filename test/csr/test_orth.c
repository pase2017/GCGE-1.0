/*
 * =====================================================================================
 *
 *       Filename:  test_cg.c
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

//#include "memwatch.h"
#include "gcge_app_csr.h"

typedef struct ORTH_WS_ {

   GCGE_DOUBLE  reorth_tol;
   GCGE_INT     max_reorth_num;
   GCGE_DOUBLE  orth_zero_tol;
   void**       vec_tmp;  /* 长度为1的多向量指针 */
   GCGE_DOUBLE* array_tmp;/* 长度为输入*end的数组指针 */

   GCGE_INT     print_level;/* 0不打印 */

}ORTH_WS;


/**
 * @brief 
 *
 * 正交化V, 从start到end, 且返回end, 因为正交化过程中会去掉线性相关的向量
 * 默认从0到start-1已经被正交化
 *
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
 * for current = (*start) ; current < (*end); ++current
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
 * @param start 从0 到 (*start)-1已经正交归一
 * @param end   输入要正交的向量个数， 返回被正交化的向量个数
 * @param B
 * @param ops
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
  */
void GCGE_Default_Orthogonalize(void **V, GCGE_INT *start, GCGE_INT *end, void *B, GCGE_OPS *ops)
{
   ORTH_WS     *workspace       = (ORTH_WS*)(ops->orthogonalize_workspace);
   GCGE_DOUBLE  reorth_tol      = workspace->reorth_tol;
   GCGE_INT     max_reorth_num  = workspace->max_reorth_num;

   GCGE_INT     print_orth_zero = workspace->print_level;
   GCGE_DOUBLE  orth_zero_tol   = workspace->orth_zero_tol;
   void       **V_tmp           = workspace->vec_tmp;
   GCGE_DOUBLE *d_tmp           = workspace->array_tmp;

   GCGE_INT     current = 0;
   GCGE_INT     reorth_count = 0;
   GCGE_INT     mv_s[2] = {0, 0};
   GCGE_INT     mv_e[2] = {0, 0};
   GCGE_DOUBLE  vin = 0.0;
   GCGE_DOUBLE  ratio = 1.0;
   void        *vec_current;
   void        *vec_tmp;
   for(current = (*start); current < (*end); ++current)
   {
      reorth_count = 0;
      do{
	 ops->GetVecFromMultiVec(V_tmp,   0, &vec_tmp);
	 ops->GetVecFromMultiVec(V, current, &vec_current);
	 if(B == NULL)
	 {
	    ops->VecAxpby(1.0, vec_current, 0.0, vec_tmp);
	 }
	 else
	 {
	    ops->MatDotVec(B, vec_current, vec_tmp);
	 }
	 ops->RestoreVecForMultiVec(V_tmp,   0, &vec_tmp);
	 ops->RestoreVecForMultiVec(V, current, &vec_current);

	 mv_s[0] = 0;         mv_s[1] = 0;
	 mv_e[0] = current+1; mv_e[1] = 1;
	 ops->MultiVecInnerProd(V, V_tmp, d_tmp, "ns", mv_s, mv_e, current+1, ops);

	 if(current != 0)
	 {
	    mv_s[0] = 0;       mv_s[1] = 0;
	    mv_e[0] = current; mv_e[1] = 1;
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

      }while(ratio < reorth_tol && reorth_count < max_reorth_num);

      if(d_tmp[current] > orth_zero_tol)
      {
	 ops->GetVecFromMultiVec(V, current, &vec_current);
	 ops->VecAxpby(0, vec_current, 1/sqrt(d_tmp[current]), vec_current);
	 ops->RestoreVecForMultiVec(V, current, &vec_current);
      }
      else 
      {
	 if(print_orth_zero == 1)
	 {
	    printf("In Orthogonal, there is a zero vector!, "
		  "current = %d, start = %d, end = %d\n", current, (*start), *end);
	 }
	 mv_s[0] = current;   mv_s[1] = *end-1;
	 mv_e[0] = current+1; mv_e[1] = *end;
	 ops->MultiVecSwap(V, V, mv_s, mv_e, ops);
	 (*end)--;
	 current--;
      }
   }
}


int main(int argc, char* argv[])
{
   //    mwInit();
   srand((unsigned)time(NULL));

   //创建矩阵
    const char *file_A = "../data/A_3.txt";
//   const char *file_A = "../data/testA";
   CSR_MAT *A = CSR_ReadMatFile(file_A);
   CSR_VEC **multi_vec, **vec_tmp;
   int num_vec = 20;
   double array_tmp[num_vec];

   GCGE_OPS *ops;
   GCGE_OPS_Create(&ops);
   GCGE_CSR_SetOps(ops);
   GCGE_OPS_Setup(ops);

   ops->BuildMultiVecByMat((void *)A, (void ***)&multi_vec, num_vec, ops);
   ops->BuildMultiVecByMat((void *)A, (void ***)&vec_tmp,   num_vec, ops);

   ops->Orthogonalize = GCGE_Default_Orthogonalize;
   ORTH_WS workspace = {1e-3, 3, 1e-8, (void **)vec_tmp, array_tmp, 0};
   ops->orthogonalize_workspace = &workspace;

   int mv_s[2], mv_e[2];

   double *xAx = calloc(num_vec*num_vec, sizeof(double));
   void  *orth_mat = (void*)A;
//   void  *orth_mat = NULL;
   for (int i = 0; i < num_vec; ++i)
   {
      mv_s[0] = num_vec-1;       
      mv_e[0] = num_vec;
      ops->MultiVecSetRandomValue((void *)multi_vec, mv_s, mv_e, ops);

      mv_s[0] = 0;      
      mv_e[0] = num_vec;
      ops->Orthogonalize((void **)multi_vec, mv_s, mv_e, orth_mat, ops);
      printf ( "%d\t%d\n", mv_s[0], mv_e[0] );
      mv_s[1] = mv_s[0];
      mv_e[1] = mv_e[0];
      ops->MatDotMultiVec(orth_mat, (void **)multi_vec, (void **)vec_tmp, mv_s, mv_e, ops);
      ops->MultiVecInnerProd((void **)multi_vec, (void **)vec_tmp, xAx, "nonsym", mv_s, mv_e, num_vec, ops );

      for (int j = 0; j < num_vec; ++j)
      {
	 for (int k = 0; k < num_vec; ++k)
	 {
//	    printf ( "%0.1e  ", xAx[k*num_vec+j] );
	 }
//	 printf ( "\n" );
      }
   }
   free(xAx);


   //释放矩阵空间
   CSR_MatFree(&A);
   ops->MultiVecDestroy((void ***)&multi_vec, num_vec, ops);
   ops->MultiVecDestroy((void ***)&vec_tmp,   num_vec, ops);

   GCGE_OPS_Free(&ops);

   //    mwTerm();
   return 0;
}
