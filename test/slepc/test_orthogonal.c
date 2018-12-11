/*
 * =====================================================================================
 *
 *       Filename:  test_solver.c
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

//#include "memwatch.h"
#include <petscsys.h>
#include <petscviewer.h>
#include <petscmat.h>

#include "gcge.h"
#include "gcge_app_slepc.h"
#include "external_slepc.h"

static char help[] = "Use GCGE-SLEPc to solve an eigensystem Ax=kBx with the matrixes loaded from files.\n";

int main(int argc, char* argv[])
{
    PetscErrorCode ierr;
    SlepcInitialize(&argc,&argv,(char*)0,help);

    //读入矩阵,其中矩阵P为求解线性方程组时的预条件矩阵
    Mat A, B;
    const char *file_A = "../data/FEM_Stiff_3.petsc.bin";
    const char *file_B = "../data/FEM_Mass_3.petsc.bin";
    SLEPC_ReadMatrixBinary(&A, file_A);
    SLEPC_ReadMatrixBinary(&B, file_B);

    //创建petsc_solver
    GCGE_SOLVER *slepc_solver;
    GCGE_SOLVER_Create(&slepc_solver);
    //设置PETSC结构的矩阵向量操作
    GCGE_SOLVER_SetSLEPCOps(slepc_solver);

    //设置一些参数
    int nev = 9;
    GCGE_SOLVER_SetNumEigen(slepc_solver, nev);//设置特征值个数

    double *eval = (double *)calloc(nev, sizeof(double)); 
    BV evec;
    slepc_solver->ops->BuildMultiVecByMat((void*)A, (void***)(&evec), nev, slepc_solver->ops);

    //将矩阵与特征对设置到petsc_solver中
    GCGE_SOLVER_SetMatA(slepc_solver, A);
    GCGE_SOLVER_SetMatB(slepc_solver, B);
    GCGE_SOLVER_SetEigenvalues(slepc_solver, eval);
    GCGE_SOLVER_SetEigenvectors(slepc_solver, (void**)evec);
    slepc_solver->workspace->evec = (void**)evec;

    //对petsc_solver进行setup，检查参数，分配工作空间等
    GCGE_SOLVER_Setup(slepc_solver);

    int x_length = 3;
    int xp_length = 2*x_length;
    int xpw_length = 3*x_length;

    slepc_solver->workspace->dim_x = x_length;
    slepc_solver->workspace->dim_xp = xp_length;
    slepc_solver->para->nev = x_length;

    void **V = slepc_solver->workspace->V;
    slepc_solver->ops->MultiVecSetRandomValue(V, 0, xpw_length, slepc_solver->ops);

    GCGE_Orthogonal(V, 0, &xp_length, (void *)B, 
	    slepc_solver->ops, slepc_solver->para, 
        slepc_solver->workspace->V_tmp, slepc_solver->workspace->subspace_dtmp);

    double *prod = (double*)calloc(nev*nev, sizeof(double));
    int i = 0;
    int j = 0;
    int mv_s[2] = {0,0};
    int mv_e[2] = {0,0};
    mv_s[0] = 0;
    mv_e[0] = xpw_length;
    mv_s[1] = 0;
    mv_e[1] = xpw_length;
    slepc_solver->ops->MatDotMultiVec((void*)B, V, (void**)evec, 
            mv_s, mv_e, slepc_solver->ops);
    slepc_solver->ops->MultiVecInnerProd(V, (void**)evec, 
            prod, "sym", mv_s, mv_e, xpw_length, slepc_solver->ops);
    printf("prod:\n");
    for(i=0; i<xpw_length; i++)
    {
        for(j=0; j<xpw_length; j++)
        {
            GCGE_Printf("%f\t", prod[i*xpw_length+j]);
        }
        printf("\n");
    }

    GCGE_BOrthogonal(V, xp_length, &xpw_length, (void *)B, 
	    slepc_solver->ops, slepc_solver->para, slepc_solver->workspace);

    slepc_solver->ops->MatDotMultiVec((void*)B, V, (void**)evec, 
            mv_s, mv_e, slepc_solver->ops);
    slepc_solver->ops->MultiVecInnerProd(V, (void**)evec, 
            prod, "sym", mv_s, mv_e, xpw_length, slepc_solver->ops);
    printf("prod:\n");
    for(i=0; i<xpw_length; i++)
    {
        for(j=0; j<xpw_length; j++)
        {
            GCGE_Printf("%f\t", prod[i*xpw_length+j]);
        }
        printf("\n");
    }

    free(prod); prod = NULL;
    //释放特征值、特征向量、KSP空间
    free(eval); eval = NULL;
    slepc_solver->ops->MultiVecDestroy((void***)(&evec), nev, slepc_solver->ops);

    //释放petsc_solver中非用户创建的空间
    GCGE_SOLVER_Free(&slepc_solver);

    //释放矩阵空间
    ierr = MatDestroy(&A);
    ierr = MatDestroy(&B);

    ierr = SlepcFinalize();

    return ierr;
}
