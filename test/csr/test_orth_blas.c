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

#define BLASname(name) d##name##_


extern double BLASname(dot) (int*, double*, int*, double*, int*);
extern double BLASname(nrm2)(int*, double*, int*);
extern void   BLASname(copy)(int*, double*, int*, double*, int*);
extern void   BLASname(axpy)(int*, double*, double*, int*, double*, int*);
extern void   BLASname(scal)(int*, double*, double*, int*);
extern void   BLASname(gemv)(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);


void Orthogonalize(double *V, int ldV, int *start, int *end)
{
   double  reorth_tol      = 1e-3;
   int     max_reorth_num  = 3;
   int     print_orth_zero = 0;
   double  orth_zero_tol   = 1e-8;

   int     current = 0;
   int     reorth_count = 0;
   double  vin = 0.0;
   double  vout = 0.0;
   double  ratio = 0.0;
   int     idx = 0;
   double  ip = 1.0;
   double  alpha= 0.0, beta = 0.0;; 

   int    inc = 1;

   /* For BLAS level 2 */
   double *Ax = calloc(*end, sizeof(double));
   /* TODO 是否可以再次利用BLAS level 3进行块的正交化 */
   for(current = (*start); current < (*end); ++current)
   {
      vout = BLASname(nrm2)(&ldV, V+current*ldV, &inc);
      if(current != 0)
      {
	 reorth_count = 0;
	 ratio = 0;
	 while(ratio < reorth_tol && reorth_count < max_reorth_num && vout > orth_zero_tol)
	 {
	    vin = vout;
	    /* BLAS level 1 */
//	    for(idx = 0; idx < current; idx++)
//	    {
//	       ip = BLASname(dot)(&ldV, V+idx*ldV, &inc, V+current*ldV, &inc);
//	       ip *= -1;
//	       BLASname(axpy)(&ldV, &ip, V+idx*ldV, &inc, V+current*ldV, &inc);
//	    }
	    /* BLAS level 2 */
	    alpha= 1.0;  beta = 0.0;
	    BLASname(gemv)("T", &ldV, &current, &alpha, V, &ldV, V+current*ldV, &inc, &beta, Ax, &inc);
	    alpha= -1.0; beta = 1.0;
	    BLASname(gemv)("N", &ldV, &current, &alpha, V, &ldV, Ax, &inc, &beta, V+current*ldV, &inc);

	    vout = BLASname(nrm2)(&ldV, V+current*ldV, &inc);
	    ratio = vout/vin;
	    reorth_count++;
	 }
      }

      if(vout > orth_zero_tol)
      {
	 vout = 1.0/vout;
	 BLASname(scal)(&ldV, &vout, V+current*ldV, &inc);
      }
      else 
      {
	 if(print_orth_zero == 1)
	 {
	    printf("In OrthogonalSubspace, there is a zero vector!, "
		  "current = %d, start = %d, end = %d\n", current, *start, *end);
	 }
	 BLASname(copy)(&ldV, V+(*end -1)*ldV, &inc, V+current*ldV, &inc);
	 (*end)--;
	 current--;
      }
   }
   free(Ax);
}


int main(int argc, char* argv[])
{
    srand((unsigned)time(NULL));
    //创建矩阵

    int n = 10;
    double entries[n*n];
    double ip[n*n];
    int ldV = n;
    int start = 0;
    int end = n;
    int inc = 1;
    int i, j, k;

    for( i=0; i<n*n; i++ )
    {
       entries[i] = 0;
    }

    start = 0;
    for( i=0; i<n; i++ )
    {
       for( j=0; j<n; j++ )
       {
	  entries[n*(n-1)+j] = (double)rand()/((double)RAND_MAX+1);
       }
       end = n;
       Orthogonalize(entries, ldV, &start, &end);
       printf ( "end = %d\n", end );
       start = end;
       printf ( "entries\n" );
       for (j = 0; j < n; ++j)
       {
	  for (k = 0; k < n; ++k)
	  {
	     printf ( "%e  ", entries[k*n+j] );
	  }
	  printf ( "\n" );
       }

       for(j=0; j<n; j++)
       {
	  for(k=0; k<n; k++)
	  {
	     ip[j*n+k] = BLASname(dot)(&n, entries+j*n, &inc, entries+k*n, &inc);
	  }
       }
       printf ( "ip\n" );
       for (j = 0; j < n; ++j)
       {
	  for (k = 0; k < n; ++k)
	  {
	     printf ( "%e  ", ip[k*n+j] );
	  }
	  printf ( "\n" );
       }
    }

    return 0;
}
