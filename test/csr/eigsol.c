/*
 * =====================================================================================
 *
 *       Filename:  gcge_eigsol.c
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

#include "gcge_eigsol.h"


typedef struct EIGSOL_WS_ {

   GCGE_INT     max_size_X;
   GCGE_INT     max_iter_num;
   GCGE_INT     *unlock;
   GCGE_DOUBLE  *subspace_matrix;
   GCGE_DOUBLE  *subspace_eval;
   GCGE_DOUBLE  *subspace_evec;
   void         **V;
   void         **RitzVec;
   GCGE_INT     given_init_evec;

   GCGE_INT     size; /* 记录指针分配的长度，可以是多个size, 
                         若分配大小可以通过max_size_X计算得到，
			 则不需要这个size
		       */
}EIGSOL_WS;

/**
 * @brief 计算 A * evec = eval * B * evec 
 *    
 *    Initialization: 
 *       对evec设定随机初值, 把evec复制给X
 *       对X正交化,在X张成的子空间中计算RR问题
 *       
 *    Loop: 
 *       对未收敛的X计算P
 *       对未收敛的X计算W
 *       对W进行正交化
 *       在[X P W]张成的子空间中计算RR问题
 *       检查特征值重数，调整X的长度
 *       检查收敛性, 得到block_size, unlock, num_unlock
 *
 *    Finalization:
 *       把X复制给evec
 *       把subspace_eval复制给eval
 *
 * @param A
 * @param B
 * @param eval
 * @param evec
 * @param nev  要求的特征值个数
 * @param ops
 */
void GCGE_Default_SolveEigenvalues(void *A, void *B, GCGE_DOUBLE *eval, void **evec, GCGE_INT nev, GCGE_OPS *ops)
{
   EIGSOL_WS   *workspace       = ops->solve_eigenvalues_workspace;

   GCGE_INT    max_size_X       = workspace->max_size_X;      /* 最大的size_X, 即最多求解Ritz值的个数  */
   GCGE_INT    max_iter_num     = workspace->max_iter_num;    /* 最大迭代次数 */

   GCGE_INT    *unlock          = workspace->unlock;          /* 存储未收敛的特征对编号, 分配的长度为? */
   GCGE_DOUBLE *subspace_matrix = workspace->subspace_matrix; /* V^T A V,                分配的长度为? */
   GCGE_DOUBLE *subspace_eval   = workspace->subspace_eval;   /* 子空间特征向量,         分配的长度为? */
   GCGE_DOUBLE *subspace_evec   = workspace->subspace_evec;   /* 子空间特征向量,         分配的长度为? */
   void        **V              = workspace->V;               /* 存储 [X P W],           分配的长度为? */
   void        **RitzVec        = workspace->RitzVec;         /* 存储Ritz向量,           分配的长度为? */

   GCGE_INT    size_X , size_P , size_W ; /* X P W 的长度, start_X == 0 */
   GCGE_INT    start_X, start_P, start_W; /* X P W 在V中的起始列号      */
   GCGE_INT    end_X  , end_P  , end_W  ; /* X P W 在V中的结束列号+1    */

   GCGE_INT    block_size = 0;            /* 分批计算特征对的个数       */
   GCGE_INT    num_iter   = 0;            /* 迭代次数                   */
   GCGE_INT    mv_s[2]    = {0, 0};
   GCGE_INT    mv_e[2]    = {0, 0};

   /* Initialization start 
    *
    * 对X做正交化
    */
   size_X  = nev;  size_P = 0;  size_W = 0;
   start_X = 0;     end_X = size_X;
   start_P = end_X; end_P = start_P+size_P;
   start_W = end_P; end_W = start_W+size_W;
   size_V  = end_W;
   /*
    * 如果用户不给初值，这里给随机初值，用户需要提供给向量设随机值的函数
    */
   if(workspace->given_init_evec == 0)
   {
      mv_s[0] = start_X; mv_e[0] = end_X;
      ops->MultiVecSetRandomValue(V, mv_s, mv_e, ops);
   }
   else
   {
      /* 把用户提供的evec初值copy给V */
      mv_s[0] = 0;       mv_e[0] = nev; 
      mv_s[1] = start_X; mv_e[1] = end_X;
      ops->MultiVecAxpby(1.0, evec, 0.0, V, mv_s, mv_e, ops);
   }
   /* 
    * 对初始近似特征向量做B正交化, 会删去线性相关的向量
    * 返回end_X
    */
   ops->Orthogonalize(V, &start_X, &end_X, B, ops);
   size_X = end_X - start_X;
   if(size_X < nev)
   {
      printf("The initial eigenvectors is linearly dependent!\n");
   }
   /* 子空间矩阵的大小 */
   size_V = end_X;

   /* Initialization end */

   /* LOOP start 
    *
    * size_V是上一次迭代[X P W]的长度, end_X是上一次的X长度
    * 基于上次迭代的[X P W]更新P, 再利用RitzVec更新X, 利用已更新的X更新W 
    * 利用新的[X P W] 更新subspace_matrix, 并求特征值和特征向量subspace_eval, subspace_evec
    * 检查重数更新size_X, 并得到RitzVec
    * 再检查收敛性更新unlock
    */
   num_unlock = nev;
   num_iter   = 0;
   while((num_unlock > 0)&&(num_iter < max_iter_num))
   {
      /*
       * 若size_V == end_X, 表明XP和XW都没有，那么不能计算P, 直接返回start_P = end_X; end_P = start_P;
       */
      if (size_V == end_X)
      {
	 start_P = end_X;
	 end_P   = start_P;
      }
      else 
      {
	 /*
	  * size_X 会被上一次迭代最后检查重数而更新, 它是Ritz向量的个数, 
	  * 此时end_X只是上一次的X个数, 即当前V中X的长度, 
	  *
	  * 计算P=X_{j}-X_{j-1}, 是通过subspace的正交操作实现P的计算以及对X的正交化
	  * subspace_evec
	  *    XX  PX                 XX  O
	  *    XP  PP  => 复制并置零  XP  XP  => 正交化
	  *    XW  PW                 XW  XW
	  *
	  * XX的行列数分别是end_X和size_X, 
	  * PX的行列数分别是end_X和block_size, 
	  *
	  * 对subspace_evec, 将unlock[j]列进行复制，并对前end_X行置零
	  * 进行subspace上的正交化之后
	  * 返回start_P = size_X, end_P = start_P + block_size ... (若没有线性相关出现)
	  * 并做线性组合, 将P更新到V的从start_P到end_P-1列
	  *
	  * 返回start_P, end_P
	  */
	 ComputeP(V, subspace_evec, end_X, size_X, &start_P, &end_P, block_size, unlock, ops);
      }
      size_P = end_P - start_P;

      /* 当第一次迭代时，只计算X生成的subspace的特征值问题 */
      start_W = end_P; 
      if (num_iter == 0)
      {
	 end_W = start_W;
      }
      else
      {
	 /*
	  * 将V中的X更新, 即V中的X部分与RitzVec交换指针
	  * 这里之所以再次赋值end_X是因为endX是上次迭代时V中的X长度
	  */
	 end_X   = firstX + size_X;
	 mv_s[0] = fisrtX; mv_e[0] = end_X;
	 mv_s[1] = 0;      mv_e[1] = size_X;
	 ops->MultiVecSwap(V, RitzVec, mv_s, mv_e, ops);
	 /*
	  * 计算AW = \lambda BX, 只计算第unlock[j]列的X对应的W, 0<=j<block_size
	  * 并将W对[X, P]做正交化, 
	  * 返回end_W, 即size_W = end_W-start_W. 这是由于过程中会将线性相关的部分去掉
	  * 若没有出现线性相关，那么size_W == block_size
	  */
	 ComputeW(V, A, B, subspace_eval, start_W, &end_W, block_size, unlock, ops);
      }
      size_W = end_W - start_W;
      /*
       * 计算子空间矩阵
       * subspace_matrix = [X, P, W]^T*A*[X, P, W]
       * subspace_matrix是size_V阶方阵，且leading dimension也是size_V
       */
      size_V = end_W;
      ComputeSubspaceMatrix(A, V, subspace_matrix, size_V, ops);
      /*
       * 求解size_X个计算子空间矩阵的特征对
       * 存储在subspace_eval, subspace_evec
       */
      ComputeSubspaceEigenpairs(subspace_matrix, size_V, subspace_eval, subspace_evec, size_X, ops);
      /*
       * 检查特征值重数, 确定size_X, 即本次更新X的个数, 
       * 完成检查重数之后, 会更改计算Ritz向量的个数, 
       * 但end_X不更新, 存储的是上一次的size_X 
       */
      CheckEvalMultiplicity(nev+2, max_size_X, &size_X, subspace_eval);
      /*
       * 用子空间特征向量与基向量V计算size_X个Ritz向量, 即
       * RitzVec[j] = [X P W] subspace_evec[j], 0<=j<size_X
       */
      ComputeX(V, subspace_evec, RitzVec, size_X, ops);
      /*
       * 检查收敛性, 计算特征对的残差，更新block_size, num_unlock和unlock
       *
       * nev是要算的特征值个数，size_X是计算出的特征值个数
       * TODO 需要加入形参tol用以检查收敛性
       * 
       * 计算
       * idx = unlock[j]  0<=j<block_size
       * A RitzVec[ idx ] - subspace_eval[ idx ] B RitzVec[ idx ] 的模并判断是否收敛
       *
       * 返回block_size, 即下次迭代时, 对unlock[j], 0<=j<block_size的X计算P和W
       * 返回unlock和num_unlock, 记录所有未收敛的特征对编号
       */
      CheckConvergence(nev, A, B, 
	    size_X, subspace_eval, RitzVec, 
	    &block_size, &num_unlock, unlock, ops);
      for(i=0; i<nev; i++)
      {
	 printf("num_iter: %d, num_unlock: %d, eval[%d]: %e\n", 
	       num_iter, num_unlock, i, subspace_eval[i]);
      }
      ++num_iter;
   }
   /* LOOP end */

   /* Finalization start 
    *
    * 把计算得到的近似特征对拷贝给eval, evec
    *
    * TODO 这里可以用BLAS进行copy特征值
    *      并且要检查是否已经收敛了nev个特征值
    *      考虑nev是否应该是一个返回值
    */
   memcpy(eval, subspace_eval, nev*sizeof(GCGE_DOUBLE));
   mv_s[0] = 0; mv_e[0] = nev;
   mv_s[1] = 0; mv_e[1] = nev;
   ops->MultiVecAxpby(1.0, RitzVec, 0.0, evec, mv_s, mv_e, ops);
   /* Finalization end */
}
