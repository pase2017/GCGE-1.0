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
    const char *file_A = "../data/A_5.txt";
//   const char *file_A = "../data/testA";
   CSR_MAT *A = CSR_ReadMatFile(file_A);
   CSR_VEC **multi_vec, **vec_tmp;
   int num_vec = 10;
   double array_tmp[num_vec];

   GCGE_OPS *ops;
   GCGE_OPS_Create(&ops);
   GCGE_CSR_SetOps(ops);
   GCGE_OPS_Setup(ops);

   ops->BuildMultiVecByMat((void *)A, (void ***)&multi_vec, num_vec, ops);
   ops->BuildMultiVecByMat((void *)A, (void ***)&vec_tmp,   num_vec, ops);

   ops->Orthogonalize = GCGE_Default_Orthogonalize;
   ORTH_WS workspace = {1e-3, 3, 1e-8, (void **)vec_tmp, array_tmp, 1};
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
	    printf ( "%e\t", xAx[k*num_vec+j] );
	 }
	 printf ( "\n" );
      }
   }
   free(xAx);


   //释放矩阵空间
   CSR_MatFree(&A);
   ops->FreeMultiVec((void ***)&multi_vec, num_vec, ops);
   ops->FreeMultiVec((void ***)&vec_tmp,   num_vec, ops);

   GCGE_OPS_Free(&ops);

   //    mwTerm();
   return 0;
}
