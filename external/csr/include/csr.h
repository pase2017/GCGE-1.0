/*
 * =====================================================================================
 *
 *       Filename:  csr.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2018年09月24日 09时50分16秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
/*
 * 定义稀疏矩阵，并且提供稀疏矩阵和相应的向量的数据结构和相互之间的操作：
 * 目前提供的操作：
 * 1. 向量的线性组合
 * 2. 向量的内积
 * 3. 矩阵向量乘
 *  创建多个向量
 */



#ifndef _CSR_H_
#define _CSR_H_

typedef struct CSR_MAT_ {	
    /** 下面是存储稀疏矩阵的行压缩形式 */
    /** number rof rows */
    int N_Rows;
    /** number columns */
    int N_Columns;
    /** number of matrix entries */
    int N_Entries;
    /** in which column is the current entry */
    int *KCol;
    /** index in KCol where each row starts */
    int *RowPtr;   
    /** matrix elements in an array */
    double *Entries;
}CSR_MAT;

typedef struct CSR_VEC_ {
    int    size;
    double *Entries;
}CSR_VEC;


//CSR矩阵向量操作的具体函数

//获取随机向量
void CSR_VecSetRandomValue(CSR_VEC *vec);
//向量的线性代数操作：y = a*x+b*y
//如果 a=0.0, 那么 y = b*y, 相当于VecScale
//如果 b=0.0, 那么 y = a*x, 相当于VecScalep2q
//如果 a=1.0, b=0.0, 那么 y = x, 相当于VecCopy
//如果 a=0.0, b=0.0, 那么 y = 0, 相当于VecSetZero
void CSR_VecAxpby(double a, CSR_VEC *x, double b, CSR_VEC *y);
//计算向量内积，如果计算向量范数，就先计算内积，再开平方
void CSR_VecInnerProd(CSR_VEC *x, CSR_VEC *y, double *xTy);
//由已给矩阵创建向量
void CSR_BuildVecByMat(CSR_MAT *mat, CSR_VEC **vec);
//由已给向量创建向量
void CSR_BuildVecByVec(CSR_VEC *s_vec, CSR_VEC **d_vec);
void CSR_BuildMultiVecByMat(CSR_MAT *mat, CSR_VEC ***vec, int nev);
//释放向量组空间
void CSR_VecFree(CSR_VEC **vec);
void CSR_MultiVecFree(CSR_VEC ***vec, int n);
//稀疏矩阵乘向量: y = mat * x
void CSR_MatDotVec(CSR_MAT *mat, CSR_VEC *x, CSR_VEC *y);
//从文件 filename 中读入稀疏矩阵
CSR_MAT *CSR_ReadMatFile(const char *filename);
void CSR_MatFree(CSR_MAT **mat);
//打印稀疏矩阵
void CSR_PrintMat(CSR_MAT *mat);
//打印向量
void CSR_PrintVec(CSR_VEC *vec);
//处理Dirichlet边条件
void CSR_SetDirichletBoundary(CSR_VEC **Vecs, int nev, CSR_MAT *A, CSR_MAT *B);
#endif
