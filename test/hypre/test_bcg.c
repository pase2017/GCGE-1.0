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

#include <math.h>
#include <time.h>
#include "HYPRE_krylov.h"

#include "gcge_app_hypre.h"
#include <float.h>

HYPRE_Int CreateLaplaceEigenProblem(MPI_Comm comm, HYPRE_Int dim, 
   HYPRE_Real** exact_eigenvalues, HYPRE_Int block_size, 
   HYPRE_Int n,
   HYPRE_IJMatrix* A, HYPRE_IJMatrix* B, 
   HYPRE_IJVector* x, HYPRE_IJVector* b );
HYPRE_Int DestroyLaplaceEigenProblem(HYPRE_Real* exact_eigenvalues, 
   HYPRE_IJMatrix A, HYPRE_IJMatrix B, 
   HYPRE_IJVector x, HYPRE_IJVector b);

int main(int argc, char* argv[])
{
   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   int myid, num_procs;
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

   /* 创建矩阵 start */
   HYPRE_IJMatrix     A,        B;
   HYPRE_IJVector     x,        b;
   HYPRE_ParCSRMatrix parcsr_A, parcsr_B;
   HYPRE_ParVector    par_x,    par_b;

//   srand((unsigned)time(NULL));
   /*----------------------- Laplace精确特征值 ---------------------*/
   /* Preliminaries: want at least one processor per row */
   int n = 33;
   int nev = 10;

   /* Parse command line */
   {
      int arg_index = 0;
      int print_usage = 0;

      while (arg_index < argc)
      {
	 if ( strcmp(argv[arg_index], "-n") == 0 )
	 {
	    arg_index++;
	    n = atoi(argv[arg_index++]);
	 }
	 else if ( strcmp(argv[arg_index], "-nev") == 0 )
	 {
	    arg_index++;
	    nev = atoi(argv[arg_index++]);
	 }
	 else if ( strcmp(argv[arg_index], "-help") == 0 )
	 {
	    print_usage = 1;
	    break;
	 }
	 else
	 {
	    arg_index++;
	 }
      }

      if ((print_usage) && (myid == 0))
      {
	 printf("\n");
	 printf("Usage: %s [<options>]\n", argv[0]);
	 printf("\n");
	 printf("  -n <n>        : problem size in each direction (default: 10)\n");
	 printf("  -nev <n>      : eigenproblem nev (default: 3)\n");
	 printf("\n");
      }

      if (print_usage)
      {
	 MPI_Finalize();
	 return (0);
      }
   }

   if (n*n < num_procs) n = sqrt(num_procs) + 1;
   int N = n*n; /* global number of rows */
   double h = 1.0/(n+1); /* mesh size*/
   /* Create 2D-laplace matrix and eigenvalues */
   double *exact_eigenvalues;
   CreateLaplaceEigenProblem(MPI_COMM_WORLD, 2, &exact_eigenvalues, nev, n, 
	 &A, &B, &x, &b);
   /* Get the parcsr matrix object to use */
   HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
   HYPRE_IJMatrixGetObject(B, (void**) &parcsr_B);
   /* Create sample rhs and solution vectors */
   HYPRE_IJVectorGetObject(x, (void **) &par_x);
   HYPRE_IJVectorGetObject(b, (void **) &par_b);
   /* 创建矩阵 end */

   /* 创建特征值求解器 */
   GCGE_SOLVER *hypre_solver = GCGE_HYPRE_Solver_Init(
	 parcsr_A, parcsr_B, nev, argc, argv);   

   /* 一些参数的设置 */
   hypre_solver->para->dirichlet_boundary = 0;
   hypre_solver->para->ev_tol             = 1e-2;
   hypre_solver->para->cg_rate            = 0.0;
   hypre_solver->para->cg_max_it          = 20;

   /* 求解特征值问题 */
//   GCGE_SOLVER_Solve(hypre_solver);  


   /* It should be modified */
   hypre_solver->workspace->evec = hypre_solver->evec;

   HYPRE_Real      *eval;
   GCGE_SOLVER_GetEigenvalues (hypre_solver, &eval);


   int seed = 314;
   int i;

   void **RHS = hypre_solver->workspace->V;
   void **X   = hypre_solver->workspace->V+nev;
   void **rhs = hypre_solver->workspace->V+2*nev;

   for (i = 0; i < nev; ++i)
   {
      HYPRE_ParVectorSetRandomValues((HYPRE_ParVector)(RHS[i]), seed);
      HYPRE_ParVectorSetRandomValues((HYPRE_ParVector)(X[i]), seed);
      seed += rand();
      HYPRE_ParVectorCopy((HYPRE_ParVector)(RHS[i]), (HYPRE_ParVector)(rhs[i]));
   }

   int global_time_index;
   global_time_index = hypre_InitializeTiming("GCGE BCG");
   hypre_BeginTiming(global_time_index);
   HYPRE_Solver precond;
   HYPRE_BoomerAMGCreate(&precond);
   /* Set some parameters */
   HYPRE_BoomerAMGSetPrintLevel(precond,  2); /* print solve info + parameters */
   HYPRE_BoomerAMGSetNumSweeps(precond,   1); /* Sweeeps on each level */
   HYPRE_BoomerAMGSetTol(precond,       0.0); /* conv. tolerance */
   HYPRE_BoomerAMGSetMaxIter(precond,    10); /* max iteration */
   HYPRE_BoomerAMGSetup(precond, parcsr_A, par_b, par_x);

   for (i = 0; i < nev; ++i)
   {
      HYPRE_BoomerAMGSolve(precond, parcsr_A, (HYPRE_ParVector)RHS[i], (HYPRE_ParVector)X[i]);
      HYPRE_ParVectorCopy((HYPRE_ParVector)(rhs[i]), (HYPRE_ParVector)(RHS[i]));
   }

   GCGE_BCG((void *)parcsr_A, RHS, X, 0, nev, 
	 hypre_solver->ops, hypre_solver->para, hypre_solver->workspace);


   hypre_EndTiming(global_time_index);
   hypre_PrintTiming("Solve phase times", MPI_COMM_WORLD);
   hypre_FinalizeTiming(global_time_index);
   hypre_ClearTiming();


   for (i = 0; i < nev; ++i)
   {
      hypre_solver->ops->VecInnerProd(X[i], X[i], eval+i);
      hypre_printf ( "%d: InnerProd = %e\n", i, eval[i] );
      hypre_solver->ops->VecLocalInnerProd(X[i], X[i], eval+i);
      hypre_printf ( "%d %d: LocalInnerProd = %e\n", i, myid, eval[i] );
   }

   double norm_rhs;
   for (i = 0; i < nev; ++i)
   {
      HYPRE_ParVectorInnerProd( (HYPRE_ParVector)rhs[i], (HYPRE_ParVector)rhs[i], &norm_rhs);
      norm_rhs = sqrt(norm_rhs);
      HYPRE_ParCSRMatrixMatvec(1.0, parcsr_A, (HYPRE_ParVector)X[i], -1.0, (HYPRE_ParVector)rhs[i] );
      HYPRE_ParVectorInnerProd( (HYPRE_ParVector)rhs[i], (HYPRE_ParVector)rhs[i], eval+i);
//      hypre_printf ( "%d: residual = %e, relative residual = %e\n", i, sqrt(eval[i]), sqrt(eval[i])/norm_rhs );
      hypre_printf ( "%d: residual = %e\n", i, sqrt(eval[i]) );
   }


   /* 释放所有空间 */
   HYPRE_BoomerAMGDestroy(precond);
   GCGE_SOLVER_Free_All(&hypre_solver);
   DestroyLaplaceEigenProblem(exact_eigenvalues, A, B, x, b);

   /* Finalize MPI*/
   MPI_Finalize();
   return 0;
}
