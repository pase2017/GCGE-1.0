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

  // srand((unsigned)time(NULL));
   srand(775);
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

   /* 创建线性解法器 start */
   int solver_id = 0;
   HYPRE_Solver linear_solver,  precond;
   if (solver_id == 0)
   {
      HYPRE_BoomerAMGCreate(&linear_solver);
      /* Set some parameters */
      HYPRE_BoomerAMGSetPrintLevel(linear_solver,  1); /* print solve info + parameters */
      HYPRE_BoomerAMGSetNumSweeps(linear_solver,   2); /* Sweeeps on each level */
      HYPRE_BoomerAMGSetTol(linear_solver,       0.0); /* conv. tolerance */
      HYPRE_BoomerAMGSetMaxIter(linear_solver,     1); /* max iteration */
      HYPRE_BoomerAMGSetup(linear_solver, parcsr_A, par_b, par_x);
   }
   else if (solver_id == 1)
   {
      /* pcg+amg */
      /*  Create solver */
      HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD,  &linear_solver);
      /*  Set some parameters */
      HYPRE_PCGSetMaxIter(linear_solver,    10); /*  max iterations */
      HYPRE_PCGSetTol(linear_solver,      1e-7); /*  conv. tolerance */
      HYPRE_PCGSetTwoNorm(linear_solver,     1); /*  use the two norm as the stopping criteria */
      HYPRE_PCGSetPrintLevel(linear_solver,  0); /*  print solve info */
      HYPRE_PCGSetLogging(linear_solver,     0); /*  needed to get run info later */
      /*  Now set up the AMG preconditioner and specify any parameters */
      HYPRE_BoomerAMGCreate(&precond);
      /* Set some parameters */
      HYPRE_BoomerAMGSetPrintLevel(precond,  1); /* print solve info + parameters */
      HYPRE_BoomerAMGSetNumSweeps(precond,   2); /* Sweeeps on each level */
      HYPRE_BoomerAMGSetTol(precond,       0.0); /* conv. tolerance */
      HYPRE_BoomerAMGSetMaxIter(precond,     1); /* max iteration */
      /*  Set the PCG preconditioner */
      HYPRE_PCGSetPrecond(linear_solver,  (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, 
	    (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,  precond);
      /*  Now setup and solve! */
      HYPRE_ParCSRPCGSetup(linear_solver, parcsr_A, par_b, par_x);
   }
   /* 创建线性解法器 end */

   /* 创建特征值求解器 solver_id: 0 AMG, 1 PCG, 2 GMRES */
   GCGE_SOLVER *hypre_solver = GCGE_HYPRE_Solver_Init_LinearSolverGivenByUser(
	 parcsr_A, parcsr_B, linear_solver, solver_id, nev, argc, argv);   
   /* 一些参数的设置 */
   hypre_solver->para->ev_max_it = 100;
   hypre_solver->para->ev_tol    = 1e-8;
   hypre_solver->para->orth_para->criterion_tol = DBL_EPSILON;
   hypre_solver->para->orth_para->orth_zero_tol = DBL_EPSILON;
/*
   hypre_solver->para->orth_para->max_reorth_time = 100;
*/
   hypre_solver->para->dirichlet_boundary = 0;
   hypre_solver->para->print_part_time = 0;
   hypre_solver->para->print_eval = 0;
   hypre_solver->para->print_para = 0;

   int global_time_index;
   global_time_index = hypre_InitializeTiming("GCGE Solver");
   hypre_BeginTiming(global_time_index);

   /* 求解特征值问题 */
   GCGE_SOLVER_Solve(hypre_solver);  

   hypre_EndTiming(global_time_index);
   hypre_PrintTiming("Solve phase times", MPI_COMM_WORLD);
   hypre_FinalizeTiming(global_time_index);
   hypre_ClearTiming();

   /* 释放所有空间 */
   GCGE_SOLVER_Free_All(&hypre_solver);
   if (solver_id == 0)
   {
      HYPRE_BoomerAMGDestroy(linear_solver);
   }
   else if (solver_id == 1)
   {
      HYPRE_ParCSRPCGDestroy(linear_solver);
      HYPRE_BoomerAMGDestroy(precond);
   }
   DestroyLaplaceEigenProblem(exact_eigenvalues, A, B, x, b);

   /* Finalize MPI*/
   MPI_Finalize();
   return 0;
}
