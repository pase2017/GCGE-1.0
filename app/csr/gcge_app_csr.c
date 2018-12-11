/*
 * =====================================================================================
 *
 *       Filename:  gcge_app_csr.c
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

/* external head file */
#include "gcge_app_csr.h"

void GCGE_CSR_SetDirichletBoundary(void**Vecs, GCGE_INT nev, void* A, void* B)
{
  CSR_SetDirichletBoundary((CSR_VEC **)Vecs, (int)nev, (CSR_MAT*)A, (CSR_MAT*)B);
}

void GCGE_CSR_VecCreateByVec(void **d_vec, void *s_vec)
{
   CSR_VecCreateByVec((CSR_VEC**)d_vec, (CSR_VEC*)s_vec);
}
void GCGE_CSR_VecCreateByMat(void **vec, void *mat)
{
   CSR_VecCreateByMat((CSR_VEC**)vec, (CSR_MAT*)mat);
}
void GCGE_CSR_VecDestroy(void **vec)
{
   CSR_VecDestroy((CSR_VEC**)vec);
}

void GCGE_CSR_VecSetRandomValue(void *vec)
{
   CSR_VecSetRandomValue((CSR_VEC*)vec);
}
void GCGE_CSR_MatDotVec(void *mat, void *x, void *r)
{
    CSR_MatDotVec((CSR_MAT*)mat, (CSR_VEC*)x, (CSR_VEC*)r);
}
void GCGE_CSR_VecAxpby(GCGE_DOUBLE a, void *x, GCGE_DOUBLE b, void *y)
{
   CSR_VecAxpby(a, (CSR_VEC*)x, b, (CSR_VEC*)y);
}
void GCGE_CSR_VecInnerProd(void *x, void *y, GCGE_DOUBLE *xTy)
{
   CSR_VecInnerProd((CSR_VEC*)x, (CSR_VEC*)y, xTy);
}
void GCGE_CSR_MultiVecPrint(void **x, GCGE_INT n)
{
   GCGE_INT i = 0;
   for(i=0; i<n; i++)
   {
    CSR_PrintVec((CSR_VEC*)x[i]); 
   }
}

void GCGE_CSR_SetOps(GCGE_OPS *ops)
{
    /* either-or */
    ops->VecCreateByVec     = GCGE_CSR_VecCreateByVec;
    ops->VecCreateByMat     = GCGE_CSR_VecCreateByMat;
    ops->VecDestroy           = GCGE_CSR_VecDestroy;

    ops->VecSetRandomValue = GCGE_CSR_VecSetRandomValue;
    ops->MatDotVec         = GCGE_CSR_MatDotVec;
    ops->VecAxpby          = GCGE_CSR_VecAxpby;
    ops->VecInnerProd      = GCGE_CSR_VecInnerProd;
    ops->VecLocalInnerProd = GCGE_CSR_VecInnerProd;
    ops->MultiVecPrint     = GCGE_CSR_MultiVecPrint;
    ops->SetDirichletBoundary = GCGE_CSR_SetDirichletBoundary;
}

void GCGE_SOLVER_SetCSROps(GCGE_SOLVER *solver)
{
    GCGE_CSR_SetOps(solver->ops);
}

//下面是一个对CSR 的GCG_Solver的初始化
GCGE_SOLVER* GCGE_CSR_Solver_Init(CSR_MAT *A, CSR_MAT *B, int num_eigenvalues, int argc, char* argv[])
{
    //第一步: 定义相应的ops,para以及workspace

    //创建一个solver变量
    GCGE_SOLVER *csr_solver;
    GCGE_SOLVER_Create(&csr_solver);
    if(num_eigenvalues != -1)
        csr_solver->para->nev = num_eigenvalues;
    //设置初始值
    GCGE_INT error = GCGE_PARA_SetFromCommandLine(csr_solver->para, argc, argv);
    int nev = csr_solver->para->nev;
    double *eval = (double *)calloc(nev, sizeof(double)); 
    csr_solver->eval = eval;
    CSR_VEC **evec = (CSR_VEC**)malloc(nev*sizeof(CSR_VEC*));
    int i = 0; 
    //这里为什么不一次性生成 CSR_MultiVecCreateByMat?
    for(i=0; i<nev; i++)
    {
        CSR_VecCreateByMat(evec+i, A);
    }

    GCGE_SOLVER_SetMatA(csr_solver, A);
    if(B != NULL)
        GCGE_SOLVER_SetMatB(csr_solver, B);
    GCGE_SOLVER_SetCSROps(csr_solver);
    GCGE_SOLVER_SetEigenvalues(csr_solver, eval);
    GCGE_SOLVER_SetEigenvectors(csr_solver, (void**)evec);
    //setup and solve
    GCGE_SOLVER_Setup(csr_solver);
    printf("Set up finish!\n");
    return csr_solver;
}





void CSR_SetDirichletBoundary(CSR_VEC **Vecs, int nev, CSR_MAT *A, CSR_MAT *B)
{
   //检查矩阵A或者矩阵B是否都是输入的矩阵
  if((A == NULL)|| (B == NULL))
     return;
  int i, j, N_Rows = A->N_Rows, start, end, length, ind_A, ind_B;
  int *A_RowPtr = A->RowPtr, *A_KCol = A->KCol, 
       *B_RowPtr = B->RowPtr, *B_KCol = B->KCol;
  double *A_Entries = A->Entries, *B_Entries = B->Entries;
  double tol = 1e-20;  

  for(i=0;i<N_Rows;i++)
  {
    //检查矩阵B的i行是否只有一个1
   start = B_RowPtr[i];
   end   = B_RowPtr[i+1];
   //ind=0是表示这一行有非零元
   ind_B = 1;
   for(j=start;j<end;j++)
   {
     if(fabs(B_Entries[j]) > tol)
     {
       ind_B = 0; 
       break;         
     }//end if 
    } //end for j
    //如果矩阵B的这一行全为零，那么去处理A的相应行
    if(ind_B == 1)
    {
      //ind_A用来记录有多少个非零元
      ind_A = 0;
      start = A_RowPtr[i];
      end   = A_RowPtr[i+1];
      for(j=start;j<end;j++)
      {
	if(fabs(A_Entries[j]) > tol)
	{
	  ind_A ++;
	}//end if 
      }//end for j
      if(ind_A == 1)
      {
	//处理边条件  
	//printf("The %d-th position is boundary!\n",i);
	for(j=0;j<nev;j++)
	{
	  //把第i个元素设置为零
	  
	  Vecs[j]->Entries[i] = 0.0;
	}//end for j=1:nev
      }//end for ind_A==1      
    }//end if ind_B ==1    
  }//end for i  
}
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
      
    int nrow, ncol, nnz, flag;
    flag = fscanf(file, "%d\n", &nrow);
    flag = fscanf(file, "%d\n", &ncol);
    flag = fscanf(file, "%d\n", &nnz);
    
    if(status) 
        printf( "(%10d, %10d, %10d)\n", nrow, ncol, nnz );
   
    int *ia = (int *)malloc( (nrow+1)*sizeof(int) );
    int *ja = (int *)malloc( nnz*sizeof(int) );
    double *aa = (double *)malloc( nnz*sizeof(double) );
    
    int i;
   
    for(i=0;i<nrow+1;i++)
    {
       flag = fscanf(file, "%d\n", ia+i);
    }
    for(i=0;i<nnz;i++)
    {
        flag = fscanf(file, "%d\n", ja+i);
    }
    for(i=0;i<nnz;i++)
    {
        flag = fscanf(file, "%lf\n", aa+i);
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
//打印稀疏矩阵
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
//打印向量
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

//随机产生向量作为初始值
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
    //把r设置位零向量。
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
void CSR_VecCreateByVec(CSR_VEC **d_vec, CSR_VEC *s_vec)
{
    int size = s_vec->size;
    (*d_vec) = (CSR_VEC *)malloc(sizeof(CSR_VEC));
    (*d_vec)->size = size;
    (*d_vec)->Entries = (double *)calloc(size, sizeof(double));
}

//由已给矩阵创建向量
void CSR_VecCreateByMat(CSR_VEC **vec, CSR_MAT *mat)
{
    int size = mat->N_Rows;
    (*vec) = (CSR_VEC *)malloc(sizeof(CSR_VEC));
    (*vec)->size = size;
    (*vec)->Entries = (double *)calloc(size, sizeof(double));
}

//由已给矩阵创建向量组
void CSR_MultiVecCreateByMat(CSR_VEC ***vec, CSR_MAT *mat, int nev)
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
void CSR_VecDestroy(CSR_VEC **vec)
{
    free((*vec)->Entries); (*vec)->Entries = NULL;
    free(*vec); (*vec) = NULL;
}

//释放向量组空间
void CSR_MultiVecDestroy(CSR_VEC ***vec, int n)
{
    int i = 0;
    for(i=0; i<n; i++)
    {
        free((*vec)[i]->Entries); (*vec)[i]->Entries = NULL;
        (*vec)[i] = NULL;
    }
    *vec = NULL;
}



