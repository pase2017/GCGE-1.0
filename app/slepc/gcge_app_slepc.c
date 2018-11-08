/*
 * =====================================================================================
 *
 *       Filename:  gcge_app_slepc.c
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
#include "gcge_app_slepc.h"

void SLEPC_ReadMatrixBinary(Mat *A, const char *filename)
{
    PetscErrorCode ierr;
    PetscViewer    viewer;
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Getting matrix...\n"); 
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer); 
    ierr = MatCreate(PETSC_COMM_WORLD, A); 
    ierr = MatSetFromOptions(*A); 
    ierr = MatLoad(*A, viewer); 
    ierr = PetscViewerDestroy(&viewer);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Getting matrix... Done\n");
}

void SLEPC_LinearSolverCreate(KSP *ksp, Mat A, Mat T)
{
    PetscErrorCode ierr;
    ierr = KSPCreate(PETSC_COMM_WORLD,ksp);
    ierr = KSPSetOperators(*ksp,A,T);
    ierr = KSPSetInitialGuessNonzero(*ksp, 1);

    ierr = KSPSetType(*ksp, KSPCG);
    //这里的rtol应取作<=ev_tol
    //PetscErrorCode  KSPSetTolerances(KSP ksp,PetscReal rtol,PetscReal abstol,PetscReal dtol,PetscInt maxits)
    ierr = KSPSetTolerances(*ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, 1000);
    PC pc;
    ierr = KSPGetPC(*ksp, &pc);
    ierr = PCSetType(pc, PCHYPRE);
    ierr = PCHYPRESetType(pc, "boomeramg");
    //最后从命令行设置参数
    ierr = KSPSetFromOptions(*ksp);
}

void SLEPC_VecLocalInnerProd(Vec x, Vec y, double *value)
{
    PetscErrorCode     ierr;
    const PetscScalar *local_x;
    const PetscScalar *local_y;
    PetscInt           low, high, length, i = 0;
    ierr = VecGetOwnershipRange(x, &low, &high);
    length = high-low;
    ierr = VecGetArrayRead(x,&local_x);
    ierr = VecGetArrayRead(y,&local_y);
    *value = 0.0;
    for(i=0; i<length; i++)
    {
        *value += local_x[i]*local_y[i];
    }
    ierr = VecRestoreArrayRead(x,&local_x);
    ierr = VecRestoreArrayRead(y,&local_y);
}


void GCGE_SLEPC_BuildVecByVec(void *s_vec, void **d_vec)
{
	PetscErrorCode ierr;
    ierr = VecDuplicate((Vec)s_vec, (Vec*)d_vec);
}
void GCGE_SLEPC_BuildVecByMat(void *mat, void **vec)
{
	PetscErrorCode ierr;
    ierr = MatCreateVecs((Mat)mat, NULL, (Vec*)vec);
}
void GCGE_SLEPC_VecFree(void **vec)
{
	PetscErrorCode ierr;
    ierr = VecDestroy((Vec*)vec);
}

void GCGE_SLEPC_VecSetRandomValue(void *vec)
{
    PetscRandom    rctx;
    PetscErrorCode ierr;
    ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
    ierr = PetscRandomSetFromOptions(rctx);
    ierr = VecSetRandom((Vec)vec, rctx);
    ierr = PetscRandomDestroy(&rctx);
}
void GCGE_SLEPC_MultiVecSetRandomValue(void **multi_vec, GCGE_INT start, GCGE_INT n_vec, struct GCGE_OPS_ *ops)
{
    PetscRandom    rctx;
    PetscErrorCode ierr;
    PetscInt       i;
    Vec            x;
    ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
    ierr = PetscRandomSetFromOptions(rctx);
    for(i=0; i<n_vec; i++)
    {
        ierr = BVGetColumn((BV)multi_vec, i, &x);
        ierr = VecSetRandom(x, rctx);
        ierr = BVRestoreColumn((BV)multi_vec, i, &x);
    }
    ierr = PetscRandomDestroy(&rctx);
}

void GCGE_SLEPC_MatDotVec(void *mat, void *x, void *r)
{
    PetscErrorCode ierr = MatMult((Mat)mat, (Vec)x, (Vec)r);
}
void GCGE_SLEPC_VecAxpby(GCGE_DOUBLE a, void *x, GCGE_DOUBLE b, void *y)
{
    PetscErrorCode ierr;
	ierr = VecScale((Vec)y, b);
    if(x != y)
    {
	    ierr = VecAXPY((Vec)y, a, (Vec)x);
    }
}
void GCGE_SLEPC_VecInnerProd(void *x, void *y, GCGE_DOUBLE *xTy)
{
	PetscErrorCode ierr = VecDot((Vec)x, (Vec)y, xTy);
}

void GCGE_SLEPC_VecLocalInnerProd(void *x, void *y, GCGE_DOUBLE *xTy)
{
    SLEPC_VecLocalInnerProd((Vec)x, (Vec)y, xTy);
}

void GCGE_SLEPC_LinearSolver(void *Matrix, void *b, void *x, struct GCGE_OPS_ *ops)
{
    PetscErrorCode ierr;
    ierr = KSPSolve((KSP)(ops->linear_solver_workspace), (Vec)b, (Vec)x);
}

void GCGE_SLEPC_GetVecFromMultiVec(void **V, GCGE_INT j, void **x)
{
    PetscErrorCode ierr = BVGetColumn((BV)V, j, (Vec*)x);
}

void GCGE_SLEPC_RestoreVecForMultiVec(void **V, GCGE_INT j, void **x)
{
    PetscErrorCode ierr = BVRestoreColumn((BV)V, j, (Vec*)x);
}

void GCGE_SLEPC_BuildMultiVecByMat(void *mat, void ***multi_vec, 
        GCGE_INT n_vec, struct GCGE_OPS_ *ops)
{
	PetscErrorCode ierr;
	Vec            vector;
    ierr = MatCreateVecs((Mat)mat, NULL, &vector);
    ierr = BVCreate(PETSC_COMM_WORLD, (BV*)multi_vec);
    ierr = BVSetType((BV)(*multi_vec), BVMAT);
    ierr = BVSetSizesFromVec((BV)(*multi_vec), vector, n_vec);
    ierr = VecDestroy(&vector);
}

void GCGE_SLEPC_FreeMultiVec(void ***MultiVec, GCGE_INT n_vec, struct GCGE_OPS_ *ops)
{
    PetscErrorCode ierr = BVDestroy((BV*)MultiVec);
}

void GCGE_SLEPC_MatDotMultiVec(void *mat, void **x, void **y, 
        GCGE_INT *start, GCGE_INT *end, struct GCGE_OPS_ *ops)
{
    PetscErrorCode ierr;
    ierr = BVSetActiveColumns((BV)x, start[0], end[0]);
    ierr = BVSetActiveColumns((BV)y, start[1], end[1]);
    ierr = BVMatMult((BV)x, (Mat)mat, (BV)y);
}

/* vec_y[j] = \sum_{i=sx}^{ex} vec_x[i] a[i][j] */
//对连续的向量组进行线性组合: y(:,start[1]:end[1]) = x(:,start[0]:end[0])*a
void GCGE_SLEPC_MultiVecLinearComb(void **x, void **y, 
        GCGE_INT *start, GCGE_INT *end,
        GCGE_DOUBLE *a, GCGE_INT lda, 
        void *dmat, GCGE_INT lddmat, struct GCGE_OPS_ *ops)
{
    PetscErrorCode ierr;
    ierr = BVSetActiveColumns((BV)x, start[0], end[0]);
    ierr = BVSetActiveColumns((BV)y, start[1], end[1]);
    Mat dense_mat;
    //y = x * dense_mat，因此dense_mat的行数为x的列数,列数为y的列数
    GCGE_INT nrows = end[0];
    GCGE_INT ncols = end[1];
    ierr = MatCreateSeqDense(PETSC_COMM_SELF, nrows, ncols, NULL, &dense_mat);
    //将稠密矩阵a中的元素赋值给dense_mat
    GCGE_DOUBLE *q;
    ierr = MatDenseGetArray(dense_mat, &q);
    memset(q, 0.0, nrows*ncols*sizeof(GCGE_DOUBLE));
    GCGE_INT i = 0;
    for(i=start[1]; i<end[1]; i++)
    {
        //默认输入的矩阵a是要从第一个元素开始用的,所以q的位置要+start[0]
        memcpy(q+i*nrows+start[0], a+(i-start[1])*lda, (end[0]-start[0])*sizeof(GCGE_DOUBLE));
    }
    ierr = MatDenseRestoreArray(dense_mat, &q);
    ierr = BVMult((BV)y, 1.0, 0.0, (BV)x, dense_mat);
    ierr = MatDestroy(&dense_mat);

}

// lda : leading dimension of matrix a
// a(start[0]:end[0],start[1]:end[1]) = V(:,start[0]:end[0])^T * W(:,start[1]:end[1])
void GCGE_SLEPC_MultiVecInnerProd(void **V, void **W, GCGE_DOUBLE *a, 
        char *is_sym, GCGE_INT *start, GCGE_INT *end, 
        GCGE_INT lda, struct GCGE_OPS_ *ops)
{
    PetscErrorCode ierr;

    ierr = BVSetActiveColumns((BV)V, start[0], end[0]);
    ierr = BVSetActiveColumns((BV)W, start[1], end[1]);
    //计算VT*W的L2内积,要先把W的矩阵设为NULL
    ierr = BVSetMatrix((BV)W, NULL, PETSC_TRUE);
    Mat dense_mat;
    //dense_mat = VT*W，因此dense_mat的行数为V的列数,列数为W的列数
    GCGE_INT nrows = end[0];
    GCGE_INT ncols = end[1];
    //col_length表示实际用到的稠密矩阵的列数,即W的列数
    GCGE_INT col_length = end[1]-start[1];
    GCGE_INT row_length = end[0]-start[0];
    ierr = MatCreateSeqDense(PETSC_COMM_SELF, nrows, ncols, NULL, &dense_mat);
    ierr = BVDot((BV)W, (BV)V, dense_mat);
    //将稠密矩阵a中的元素赋值给dense_mat
    GCGE_DOUBLE *q;
    ierr = MatDenseGetArray(dense_mat, &q);
    GCGE_INT i = 0;

    for(i=0; i<col_length; i++)
    {
        //默认输入的矩阵a是要从第一个元素开始用的,q的位置要从start[1]列开始用,每列加start[0]
        memcpy(a+i*lda, q+(start[1]+i)*nrows+start[0], row_length*sizeof(GCGE_DOUBLE));
    }
    ierr = MatDenseRestoreArray(dense_mat, &q);
    ierr = MatDestroy(&dense_mat);

}

//计算 y = a * x + b * y
void GCGE_SLEPC_MultiVecAxpby(GCGE_DOUBLE a, void **x, GCGE_DOUBLE 
        b, void **y, GCGE_INT *start, GCGE_INT *end, struct GCGE_OPS_ *ops)
{
    PetscErrorCode ierr;
    ierr = BVSetActiveColumns((BV)x, start[0], end[0]);
    ierr = BVSetActiveColumns((BV)y, start[1], end[1]);
    //If matrix Q is NULL, then an AXPY operation Y = beta*Y + alpha*X is done
    ierr = BVMult((BV)y, a, b, (BV)x, NULL);
}

//把V_1和V_2的相应的列组交换: 即： V_1(:,start[0]:end[0])与V_2(:,start[1]:end[1])相互交换
void GCGE_SLEPC_MultiVecSwap(void **V_1, void **V_2, 
        GCGE_INT *start, GCGE_INT *end, struct GCGE_OPS_ *ops)
{
    PetscErrorCode ierr;
    if(V_1 != V_2)
    {
        ierr = BVSetActiveColumns((BV)V_1, start[0], end[0]);
        ierr = BVSetActiveColumns((BV)V_2, start[1], end[1]);
        //BVCopy是将前面的copy给后面的
        ierr = BVCopy((BV)V_2, (BV)V_1);
    }
    else
    {
        Vec x, y; 
        GCGE_INT i = 0;
        GCGE_INT length = end[0]-start[0];

        for(i=0; i<length; i++)
        {
            ierr = BVGetColumn((BV)V_1, i+start[0], &x);
            ierr = BVGetColumn((BV)V_2, i+start[1], &y);
            //VecCopy是将前面的copy给后面的
            ierr = VecCopy(y, x);
            ierr = BVRestoreColumn((BV)V_1, i+start[0], &x);
            ierr = BVRestoreColumn((BV)V_2, i+start[1], &y);
        }
    }
}

void GCGE_SLEPC_SetOps(GCGE_OPS *ops)
{
    /* either-or */
    ops->BuildVecByVec     = GCGE_SLEPC_BuildVecByVec;
    ops->BuildVecByMat     = GCGE_SLEPC_BuildVecByMat;
    ops->FreeVec           = GCGE_SLEPC_VecFree;

    ops->VecSetRandomValue = GCGE_SLEPC_VecSetRandomValue;
    ops->MatDotVec         = GCGE_SLEPC_MatDotVec;
    ops->VecAxpby          = GCGE_SLEPC_VecAxpby;
    ops->VecInnerProd      = GCGE_SLEPC_VecInnerProd;
    ops->VecLocalInnerProd = GCGE_SLEPC_VecLocalInnerProd;

    ops->FreeMultiVec           = GCGE_SLEPC_FreeMultiVec;
    ops->BuildMultiVecByMat     = GCGE_SLEPC_BuildMultiVecByMat;
    ops->MultiVecSetRandomValue = GCGE_SLEPC_MultiVecSetRandomValue;

    ops->GetVecFromMultiVec = GCGE_SLEPC_GetVecFromMultiVec;
    ops->RestoreVecForMultiVec = GCGE_SLEPC_RestoreVecForMultiVec;
    ops->MatDotMultiVec = GCGE_SLEPC_MatDotMultiVec;

    ops->MultiVecLinearComb = GCGE_SLEPC_MultiVecLinearComb;
    ops->MultiVecInnerProd = GCGE_SLEPC_MultiVecInnerProd;
    ops->MultiVecSwap = GCGE_SLEPC_MultiVecSwap;
    ops->MultiVecAxpby = GCGE_SLEPC_MultiVecAxpby;

}

void GCGE_SOLVER_SetSLEPCOps(GCGE_SOLVER *solver)
{
    GCGE_SLEPC_SetOps(solver->ops);
}

//下面是一个对SLEPC 的GCG_Solver的初始化
//用户只提供要求解特征值问题的矩阵A,B,线性解法器等都用默认值
//如果用户想在命令行提供特征值个数，需要将num_eigenvalues赋值为-1
GCGE_SOLVER* GCGE_SLEPC_Solver_Init(Mat A, Mat B, int num_eigenvalues, int argc, char* argv[])
{
    //第一步: 定义相应的ops,para以及workspace

    //创建一个solver变量
    GCGE_SOLVER *slepc_solver;
    GCGE_SOLVER_Create(&slepc_solver);
    if(num_eigenvalues != -1)
        slepc_solver->para->nev = num_eigenvalues;
    GCGE_INT error = GCGE_PARA_SetFromCommandLine(slepc_solver->para, argc, argv);
    //设置初始值
    int nev = slepc_solver->para->nev;
    double *eval = (double *)calloc(nev, sizeof(double)); 
    slepc_solver->eval = eval;

	PetscErrorCode ierr;
    Vec *evec;
	Vec vector;
    ierr = MatCreateVecs(A, NULL, &vector);
    ierr = VecDuplicateVecs(vector, nev, &evec);

    GCGE_SOLVER_SetMatA(slepc_solver, A);
    if(B != NULL)
        GCGE_SOLVER_SetMatB(slepc_solver, B);
    GCGE_SOLVER_SetSLEPCOps(slepc_solver);
    GCGE_SOLVER_SetEigenvalues(slepc_solver, eval);
    GCGE_SOLVER_SetEigenvectors(slepc_solver, (void**)evec);
    //setup and solve
    GCGE_SOLVER_Setup(slepc_solver);
    GCGE_PrintParaInfo(slepc_solver->para);
    GCGE_Printf("Set up finish!\n");
    return slepc_solver;
}

//下面是一个对SLEPC 的GCG_Solver的初始化
//用户只提供要求解特征值问题的矩阵A,B,线性解法器使用KSPSolve,用户不接触ksp,设置参数均从命令行设置
GCGE_SOLVER* GCGE_SLEPC_Solver_Init_KSPDefault(Mat A, Mat B, Mat P, int num_eigenvalues, int argc, char* argv[])
{
    //首先用矩阵A,B创建slepc_solver
    GCGE_SOLVER *slepc_solver = GCGE_SLEPC_Solver_Init(A, B, num_eigenvalues, argc, argv);
    KSP ksp;
    //使用矩阵A和预条件矩阵P创建ksp,ksp的参数均为默认参数,要修改ksp参数只能用命令行参数进行设置
    SLEPC_LinearSolverCreate(&ksp, (Mat)A, (Mat)P);
    //把ksp作为ops->linear_solver_workspace
    GCGE_SOLVER_SetOpsLinearSolverWorkspace(slepc_solver, (void*)ksp);
    //将线性解法器设为KSPSolve
    slepc_solver->ops->LinearSolver = GCGE_SLEPC_LinearSolver;
    return slepc_solver;
}

//下面是一个对SLEPC 的GCG_Solver的初始化
//用户提供要求解特征值问题的矩阵A,B,以及线性解法器结构ksp,使用KSPSolve作为线性解法器
GCGE_SOLVER* GCGE_SLEPC_Solver_Init_KSPGivenByUser(Mat A, Mat B, KSP ksp, int num_eigenvalues, int argc, char* argv[])
{
    //首先用矩阵A,B创建slepc_solver
    GCGE_SOLVER *slepc_solver = GCGE_SLEPC_Solver_Init(A, B, num_eigenvalues, argc, argv);
    //用户已创建好ksp,直接赋给slepc_solver
    GCGE_SOLVER_SetOpsLinearSolverWorkspace(slepc_solver, (void*)ksp);
    //将线性解法器设为KSPSolve
    slepc_solver->ops->LinearSolver = GCGE_SLEPC_LinearSolver;
    return slepc_solver;
}

void GCGE_SOLVER_SetSLEPCOpsLinearSolver(GCGE_SOLVER *solver, KSP ksp)
{
    GCGE_SOLVER_SetOpsLinearSolverWorkspace(solver, (void*)ksp);
    //将线性解法器设为KSPSolve
    solver->ops->LinearSolver = GCGE_SLEPC_LinearSolver;
}
