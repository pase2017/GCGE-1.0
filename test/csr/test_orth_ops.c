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
 *     mv_s = {0, 0}
 *     mv_e = {current, 1}
 *     MultiVecLinearComb (V, &V_tmp, mv_s, mv_e, d_tmp, current, ops)
 *     VecAxpby(-1, V_tmp, 1, V[current]);
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
      void *B, GCGE_OPS *ops, GCGE_ORTH_PARA *orth_para, void *V_tmp, GCGE_DOUBLE *d_tmp)
{

   {
      //从V中start位置开始进行正交化
      //B正交化，先计算前start个向量的矩阵乘向量放到V_tmp中
      
    void (*MatDotMultiVec)          (void *mat, void **x, void **y, GCGE_INT *start, GCGE_INT *end, 
                                     struct GCGE_OPS_ *ops);

      ops->MatDotMultiVec(B, )
      for( i=0; i<start; i++ )
      {
	 ops->MatDotVec(B, V[i], V_tmp[i]);
      }
      for( i=start; i<(*end); i++ )
      {
	 //如果i==0,说明start=0,对第一个向量直接进行单位化
	 if(i == 0)
	 {
	    dd = sqrt(GCGE_VecMatrixVec(V[0], B, V[0], V_tmp[0], ops));
	    if(dd > orth_zero_tol)
	    {
	       ops->VecAxpby(0.0, V[0], 1.0/dd, V[0]);
	       Ind[n_nonzero++] = 0;
	       ops->MatDotVec(B, V[0], V_tmp[0]);
	    }
	 }
	 else
	 {
	    //vout = sqrt(GCGE_VecMatrixVec(V[i], B, V[i], V_tmp[start+n_nonzero], ops));
	    reorth_time = 0;
	    //进行最多max_reorth_time次重正交化
	    do{
	       ops->MatDotVec(B, V[i], V_tmp[start+n_nonzero]);
	       reorth_time += 1;
	       //先计算V[i]与前面的V[j]的局部内积，存放在inner_product中
	       for( j=0; j<start+n_nonzero+1; j++ )
	       {
		  ops->VecLocalInnerProd(V[i], V_tmp[j], inner_product+j);
	       }
	       //#if USE_MPI
#if 0
	       //共有start+n_nonzero+1个局部内积要发送
	       MPI_Allreduce(MPI_IN_PLACE, inner_product, start+n_nonzero+1, 
		     MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
	       //现在inner_product中为V[i]与相应的V[j]的内积，下面减去相应的分量
	       for( j=0; j<start; j++ )
	       {
		  ops->VecAxpby(-inner_product[j], V[j], 1.0, V[i]);
	       }
	       for( j=0; j<n_nonzero; j++ )
	       {
		  ops->VecAxpby(-inner_product[start+j], V[Ind[j]], 1.0, V[i]);
	       }
	       //下面对V[i]进行归一化，首先计算V[i]的范数
	       vout = inner_product[start+n_nonzero];
	       vin = sqrt(vout);
	       for( j=0; j<start+n_nonzero; j++ )
	       {
		  vout -= inner_product[j]*inner_product[j];
	       }
	       vout = sqrt(vout);
	    }while((vout/vin < reorth_tol)&&(reorth_time < max_reorth_time));
	    //如果向量非0，进行单位化
	    if(vout > orth_zero_tol)
	    {
	       ops->VecAxpby(0.0, V[i], 1.0/vout, V[i]);
	       ops->MatDotVec(B, V[i], V_tmp[start+n_nonzero]);
	       Ind[n_nonzero++] = i;
	    }
	    else
	    {
	       if(print_orthzero == 1)
	       {
		  printf("In Orthogonal, there is a zero vector!, "
			"i = %d, start = %d, end = %d\n", i, start, *end);
	       }
	       //Nonzero_Vec[n_zero++] = V[i];
	    }
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

    int num_vec = 2;
    CSR_VEC **multi_vec = (CSR_VEC**)malloc(num_vec*sizeof(CSR_VEC*));
    int i = 0; 
    for(i=0; i<num_vec; i++)
    {
        CSR_BuildVecByMat(A, multi_vec+i);
    }
    CSR_VEC *vec;
    CSR_BuildVecByMat(A, &vec);

    for(i=0; i<vec->size; i++)
    {
       vec->Entries[i] = i+1;
    }

    int k;
    for(k=0; k<num_vec; k++)
    {
       for(i=0; i<multi_vec[k]->size; i++)
       {
	  multi_vec[k]->Entries[i] = i+k;
       }
       CSR_PrintVec(multi_vec[k]);
    }
    int start = 0, end = 2;
    /* SHOULD be modified to ops */
    GCGE_Orthogonal((void **)multi_vec, start, &end, (void *)B, 
	  solver->ops, solver->para->orth_para, solver->workspace);

    for(k=0; k<num_vec; k++)
    {
       CSR_PrintVec(multi_vec[k]);
    }

    double vTv;
    solver->ops->VecInnerProd((void *)multi_vec[0], (void *)multi_vec[1], &vTv);
    printf ( "vTv = %f\n", vTv );

    double mvTmv[4];
    int    mv_s[2] = {0, 0};
    int    mv_e[2] = {2, 2};
    solver->ops->MultiVecInnerProd((void **)multi_vec, (void **)multi_vec, (double *)&mvTmv, 
	  "S", mv_s, mv_e, 2, solver->ops);

    printf ( "mvTmv\n" );
    printf ( "%f\t%f\n%f\t%f\n", mvTmv[0], mvTmv[1], mvTmv[2], mvTmv[3] );
    
    double entries[6] = {0.0, 1.0, 2.0, 1.0, 2.0, 3.0};
    end = 2;
    GCGE_OrthogonalSubspace((double *)entries, 0, &end, 3, 
	  NULL, -1, solver->para->orth_para, solver->workspace);
    printf ( "entries\n" );
    printf ( "%f\t%f\t%f\n%f\t%f\t%f\n", 
	  entries[0], entries[1], entries[2], 
	  entries[3], entries[4], entries[5]  );


    //释放矩阵空间
    GCGE_SOLVER_Free(&solver);
    CSR_MatFree(&A);
    CSR_MatFree(&B);
    for(i=0; i<num_vec; i++)
    {
        CSR_VecFree(multi_vec+i);
    }
    free(multi_vec); multi_vec = NULL;
    CSR_VecFree(&vec);

    return 0;
}
