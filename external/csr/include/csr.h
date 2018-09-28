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
//稀疏矩阵乘向量
void CSR_MatDotVec(CSR_MAT *mat, CSR_VEC *x, CSR_VEC *y);
CSR_MAT *CSR_ReadMatFile(const char *filename);
void CSR_MatFree(CSR_MAT **mat);

void CSR_PrintMat(CSR_MAT *mat);
void CSR_PrintVec(CSR_VEC *vec);

#endif
