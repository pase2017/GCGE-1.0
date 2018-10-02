/*
 * =====================================================================================
 *
 *       Filename:  csr.c
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

#include "csr.h"

//从文件读入CSR格式矩阵
CSR_MAT *CSR_ReadMatFile(const char *filename)
{
    int status = 1;
    FILE *file;
    file = fopen(filename,"r");
    if(!file)
    {
        printf("\ncannot open %s!\n", filename);
		exit(0);
    }
    if(status) printf("Read matrix: %s\n", filename); 
      
    int nrow, ncol, nnz;
    fscanf(file, "%d\n", &nrow);
    fscanf(file, "%d\n", &ncol);
    fscanf(file, "%d\n", &nnz);
    
    if(status) 
        printf( "(%10d, %10d, %10d)\n", nrow, ncol, nnz );
   
    int *ia = (int *)malloc( (nrow+1)*sizeof(int) );
    int *ja = (int *)malloc( nnz*sizeof(int) );
    double *aa = (double *)malloc( nnz*sizeof(double) );
    
    int i;
   
    for(i=0;i<nrow+1;i++)
    {
        fscanf(file, "%d\n", ia+i);
    }
    for(i=0;i<nnz;i++)
    {
        fscanf(file, "%d\n", ja+i);
    }
    for(i=0;i<nnz;i++)
    {
        fscanf(file, "%lf\n", aa+i);
    }
    
    CSR_MAT *matrix = malloc( sizeof(CSR_MAT) );
    
    matrix->N_Rows = nrow;
    matrix->N_Columns = ncol;
    matrix->N_Entries = nnz;
    
    matrix->RowPtr = ia;
    matrix->KCol = ja;
    matrix->Entries = aa;

    fclose(file);
    return matrix;
}

void CSR_PrintMat(CSR_MAT *mat)
{
   int nrow, ncol, idx;
   printf ( "num rows = %d\n", mat->N_Rows );
   printf ( "num cols = %d\n", mat->N_Columns );
   idx = 0;
   for (nrow  = 0; nrow < mat->N_Rows; ++nrow)
   {
      for (ncol  = 0; ncol < mat->N_Columns; ++ncol)
      {
	 if (mat->RowPtr[nrow] <= idx && mat->RowPtr[nrow+1] > idx &&  mat->KCol[idx] == ncol)
	 {
	    printf ( "%0.4f\t", mat->Entries[idx] );
	    ++idx;
	 }
	 else 
	 {
	    printf ( "%0.4f\t", 0.0 );
	 }
      }
      printf ( "\n" );

   }
}
void CSR_PrintVec(CSR_VEC *vec)
{
   int nrow;
   printf ( "num size = %d\n", vec->size );
   for (nrow  = 0; nrow < vec->size; ++nrow)
   {
      printf ( "%0.4f\t", vec->Entries[nrow] );
   }
   printf ( "\n" );
}

//释放矩阵空间
void CSR_MatFree(CSR_MAT **mat)
{
	free((*mat)->RowPtr);    (*mat)->RowPtr  = NULL;
	free((*mat)->KCol);      (*mat)->KCol    = NULL;
	free((*mat)->Entries);   (*mat)->Entries = NULL;
	free((*mat));            (*mat)          = NULL;
}

void CSR_VecSetRandomValue(CSR_VEC *vec)
{
    int j, size = vec->size;
    //time_t ts;

    //srand((unsigned)time(&ts));
    for( j=0; j<size; j++ )
    {
       vec->Entries[j] = ((double)rand())/((double)RAND_MAX+1);
    }
}

//稀疏矩阵乘向量
void CSR_MatDotVec(CSR_MAT *mat, CSR_VEC *x, CSR_VEC *r)
{
    if (mat == NULL)
    {
       CSR_VecAxpby(1.0, x, 0.0, r);
       return;
    }
    int    N_Rows = mat->N_Rows, 
           N_Columns = mat->N_Columns,
	   i, j, length, start, end,
	  *RowPtr,*KCol;
    double *Mat_Entries, *x_Entries, *r_Entries, tmp;
    RowPtr      = mat->RowPtr;
    KCol        = mat->KCol;
    Mat_Entries = mat->Entries;
    x_Entries   = x->Entries;
    r_Entries   = r->Entries;

    memset(r_Entries,0.0,N_Rows*sizeof(double));

    start = RowPtr[0];
    for(i=0;i<N_Rows;i++)
    {    
        end = RowPtr[i+1];
        length = end - start;
        for(j=0;j<length;j++)
        {
            r_Entries[i] += Mat_Entries[start+j]*x_Entries[KCol[start+j]];      
        }
        start = end;
    }
}

//向量的线性代数操作：y = a*x+b*y
//如果 a=0.0, 那么 y = b*y, 相当于VecScale
//如果 b=0.0, 那么 y = a*x, 相当于VecScalep2q
//如果 a=1.0, b=0.0, 那么 y = x, 相当于VecCopy
//如果 a=0.0, b=0.0, 那么 y = 0, 相当于VecSetZero
void CSR_VecAxpby(double a, CSR_VEC *x, double b, CSR_VEC *y)
{
    int    i, size = x->size;
    double *x_Entries = x->Entries, *y_Entries = y->Entries;
    if(a == 0.0)
    {
        for(i=0; i<size; i++)
        {
            y_Entries[i] = b * y_Entries[i];
        }
    }
    else if(b == 0.0)
    {
        if(a == 1.0)
        {
            for(i=0; i<size; i++)
            {
                y_Entries[i] = x_Entries[i];
            }
        }
        else
        {
            for(i=0; i<size; i++)
            {
                y_Entries[i] = a * x_Entries[i];
            }
        }
    }
    else
    {
        for(i=0; i<size; i++)
        {
            y_Entries[i] = a * x_Entries[i] + b * y_Entries[i];
        }
    }
}

//计算向量内积，如果计算向量范数，就先计算内积，再开平方
void CSR_VecInnerProd(CSR_VEC *x, CSR_VEC *y, double *xTy)
{
    int    i, size = x->size;
    double *x_Entries = x->Entries, *y_Entries = y->Entries, xy = 0.0;
    for(i=0; i<size; i++)
    {
        xy += x_Entries[i]*y_Entries[i];
    }
    //return xTy;
    *xTy = xy;
}

//由已给向量创建向量组
void CSR_BuildVecByVec(CSR_VEC *s_vec, CSR_VEC **d_vec)
{
    int size = s_vec->size;
    (*d_vec) = (CSR_VEC *)malloc(sizeof(CSR_VEC));
    (*d_vec)->size = size;
    (*d_vec)->Entries = (double *)calloc(size, sizeof(double));
}

//由已给矩阵创建向量
void CSR_BuildVecByMat(CSR_MAT *mat, CSR_VEC **vec)
{
    int size = mat->N_Rows;
    (*vec) = (CSR_VEC *)malloc(sizeof(CSR_VEC));
    (*vec)->size = size;
    (*vec)->Entries = (double *)calloc(size, sizeof(double));
}

//由已给矩阵创建向量组
void CSR_BuildMultiVecByMat(CSR_MAT *mat, CSR_VEC ***vec, int nev)
{
    int size = mat->N_Rows;
    int i = 0;
    *vec = (CSR_VEC**)malloc(nev*sizeof(CSR_VEC*));
    for(i=0; i<nev; i++)
    {
        (*vec)[i] = (CSR_VEC *)malloc(sizeof(CSR_VEC));
        (*vec)[i]->size = size;
        (*vec)[i]->Entries = (double *)calloc(size, sizeof(double));
    }
}

//释放向量组空间
void CSR_VecFree(CSR_VEC **vec)
{
    free((*vec)->Entries); (*vec)->Entries = NULL;
    free(*vec); (*vec) = NULL;
}

//释放向量组空间
void CSR_MultiVecFree(CSR_VEC ***vec, int n)
{
    int i = 0;
    for(i=0; i<n; i++)
    {
        free((*vec)[i]->Entries); (*vec)[i]->Entries = NULL;
        (*vec)[i] = NULL;
    }
    *vec = NULL;
}
