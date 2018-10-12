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

#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"


#include "gcge.h"
#include "gcge_app_hypre.h"

static int cmp( const void *a ,  const void *b )
{   return *(double *)a > *(double *)b ? 1 : -1; }

int main(int argc, char* argv[])
{
   //mwInit();
   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   int myid, num_procs;
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

   /* 创建矩阵 start */
   int i, k;
   int ilower, iupper;
   int local_size, extra;

   HYPRE_IJMatrix A;
   HYPRE_ParCSRMatrix parcsr_A;
   HYPRE_IJMatrix B;
   HYPRE_ParCSRMatrix parcsr_B;

   srand((unsigned)time(NULL));
   /*----------------------- Laplace精确特征值 ---------------------*/
   /* Preliminaries: want at least one processor per row */
   int n = 500;
   int nev = 10;
   if (n*n < num_procs) n = sqrt(num_procs) + 1;
   int N = n*n; /* global number of rows */
   double h = 1.0/(n+1); /* mesh size*/
   double h2 = h*h;

   HYPRE_Real *exact_eigenvalues;
   {
      int tmp_nn = (int) sqrt(nev) + 3;
      exact_eigenvalues = hypre_CTAlloc (HYPRE_Real, tmp_nn*tmp_nn);
      for (i = 0; i < tmp_nn; ++i) 
      {
	 for (k = 0; k < tmp_nn; ++k) 
	 {
	    //	    exact_eigenvalues[i*tmp_nn+k] = M_PI*M_PI*(pow(i+1, 2)+pow(k+1, 2));
	    exact_eigenvalues[i*tmp_nn+k] = ( 4*sin( (i+1)*M_PI/(2*(n+1)) )*sin( (i+1)*M_PI/(2*(n+1)) ) 
		  + 4*sin( (k+1)*M_PI/(2*(n+1)) )*sin( (k+1)*M_PI/(2*(n+1)) ) );
	 }
      }
      qsort(exact_eigenvalues, tmp_nn*tmp_nn, sizeof(double), cmp);
   }

   /* Each processor knows only of its own rows - the range is denoted by ilower
      and iupper.  Here we partition the rows. We account for the fact that
      N may not divide evenly by the number of processors. */
   local_size = N/num_procs;
   extra = N - local_size*num_procs;

   ilower = local_size*myid;
   ilower += hypre_min(myid, extra);

   iupper = local_size*(myid+1);
   iupper += hypre_min(myid+1, extra);
   iupper = iupper - 1;

   /* How many rows do I have? */
   local_size = iupper - ilower + 1;

   /* Create the matrix.
      Note that this is a square matrix, so we indicate the row partition
      size twice (since number of rows = number of cols) */
   HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
   HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &B);

   /* Choose a parallel csr format storage (see the User's Manual) */
   HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
   HYPRE_IJMatrixSetObjectType(B, HYPRE_PARCSR);

   /* Initialize before setting coefficients */
   HYPRE_IJMatrixInitialize(A);
   HYPRE_IJMatrixInitialize(B);

   /* Now go through my local rows and set the matrix entries.
      Each row has at most 5 entries. For example, if n=3:

      A = [M -I 0; -I M -I; 0 -I M]
      M = [4 -1 0; -1 4 -1; 0 -1 4]

      Note that here we are setting one row at a time, though
      one could set all the rows together (see the User's Manual).
   */
   {
      int nnz;
      double values[5];
      int cols[5];

      for (i = ilower; i <= iupper; i++)
      {
	 nnz = 0;

	 /* The left identity block:position i-n */
	 if ((i-n)>=0)
	 {
	    cols[nnz] = i-n;
	    values[nnz] = -1.0;
	    nnz++;
	 }

	 /* The left -1: position i-1 */
	 if (i%n)
	 {
	    cols[nnz] = i-1;
	    values[nnz] = -1.0;
	    nnz++;
	 }

	 /* Set the diagonal: position i */
	 cols[nnz] = i;
	 values[nnz] = 4.0;
	 nnz++;

	 /* The right -1: position i+1 */
	 if ((i+1)%n)
	 {
	    cols[nnz] = i+1;
	    values[nnz] = -1.0;
	    nnz++;
	 }

	 /* The right identity block:position i+n */
	 if ((i+n)< N)
	 {
	    cols[nnz] = i+n;
	    values[nnz] = -1.0;
	    nnz++;
	 }

	 /* Set the values for row i */
	 HYPRE_IJMatrixSetValues(A, 1, &nnz, &i, cols, values);
      }
   }
   {
      int nnz;
      double values[5];
      int cols[5];
      for (i = ilower; i <= iupper; i++)
      {
	 nnz = 1;
	 cols[0] = i;
	 values[0] = 1.0*h2;
	 /* Set the values for row i */
	 HYPRE_IJMatrixSetValues(B, 1, &nnz, &i, cols, values);
      }
   }
   /* Assemble after setting the coefficients */
   HYPRE_IJMatrixAssemble(A);
   HYPRE_IJMatrixAssemble(B);
   /* Get the parcsr matrix object to use */
   HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
   HYPRE_IJMatrixGetObject(B, (void**) &parcsr_B);
   /* Create par vector for linear solver */
   HYPRE_IJVector x;
   HYPRE_ParVector par_x;
   HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x);
   HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(x);
   HYPRE_IJVectorAssemble(x);
   HYPRE_IJVectorGetObject(x, (void **) &par_x);
   /* 创建矩阵 end */

   int solver_id = 1;
   HYPRE_Solver linear_solver,  precond;
   if (solver_id == 0)
   {
      /* 创建线性解法器amg */
      /* Create solver */
      HYPRE_BoomerAMGCreate(&linear_solver);
      /* Set some parameters (See Reference Manual for more parameters) */
      HYPRE_BoomerAMGSetPrintLevel(linear_solver, 0);  /* print solve info + parameters */
      HYPRE_BoomerAMGSetOldDefault(linear_solver); /* Falgout coarsening with modified classical interpolaiton */
      HYPRE_BoomerAMGSetRelaxType(linear_solver, 3);   /* G-S/Jacobi hybrid relaxation */
      HYPRE_BoomerAMGSetRelaxOrder(linear_solver, 1);   /* uses C/F relaxation */
      HYPRE_BoomerAMGSetNumSweeps(linear_solver, 1);   /* Sweeeps on each level */
      HYPRE_BoomerAMGSetMaxLevels(linear_solver, 20);  /* maximum number of levels */
      HYPRE_BoomerAMGSetTol(linear_solver, 1e-7);      /* conv. tolerance */
      /* Now setup and solve! */
      HYPRE_BoomerAMGSetup(linear_solver, parcsr_A, par_x, par_x);
   }
   else if (solver_id == 1)
   {
      /* 创建线性解法器pcg+amg */
      /*  Create solver */
      HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD,  &linear_solver);
      /*  Set some parameters (See Reference Manual for more parameters) */
      HYPRE_PCGSetMaxIter(linear_solver,  1000); /*  max iterations */
      HYPRE_PCGSetTol(linear_solver,  1e-7); /*  conv. tolerance */
      HYPRE_PCGSetTwoNorm(linear_solver,  1); /*  use the two norm as the stopping criteria */
      HYPRE_PCGSetPrintLevel(linear_solver,  0); /*  print solve info */
      HYPRE_PCGSetLogging(linear_solver,  0); /*  needed to get run info later */
      /*  Now set up the AMG preconditioner and specify any parameters */
      HYPRE_BoomerAMGCreate(&precond);
      HYPRE_BoomerAMGSetPrintLevel(precond,  0); /*  print amg solution info */
      HYPRE_BoomerAMGSetCoarsenType(precond,  6);
      HYPRE_BoomerAMGSetOldDefault(precond);
      HYPRE_BoomerAMGSetRelaxType(precond,  6); /*  Sym G.S./Jacobi hybrid */
      HYPRE_BoomerAMGSetNumSweeps(precond,  1);
      HYPRE_BoomerAMGSetTol(precond,  0.0); /*  conv. tolerance zero */
      HYPRE_BoomerAMGSetMaxIter(precond,  1); /*  do only one iteration! */
      /*  Set the PCG preconditioner */
      HYPRE_PCGSetPrecond(linear_solver,  (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, 
	    (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,  precond);
      /*  Now setup and solve! */
      HYPRE_ParCSRPCGSetup(linear_solver, parcsr_A, par_x, par_x);
   }
   else if (solver_id == 2)
   {
      /* 创建线性解法器gmres+amg */
      /* Create solver */
      HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, &linear_solver);
      /* Set some parameters (See Reference Manual for more parameters) */
      HYPRE_FlexGMRESSetKDim(linear_solver, 30);
      HYPRE_FlexGMRESSetMaxIter(linear_solver, 1000); /* max iterations */
      HYPRE_FlexGMRESSetTol(linear_solver, 1e-7); /* conv. tolerance */
      HYPRE_FlexGMRESSetPrintLevel(linear_solver, 2); /* print solve info */
      HYPRE_FlexGMRESSetLogging(linear_solver, 0); /* needed to get run info later */
      /* Now set up the AMG preconditioner and specify any parameters */
      HYPRE_BoomerAMGCreate(&precond);
      HYPRE_BoomerAMGSetPrintLevel(precond, 1); /* print amg solution info */
      HYPRE_BoomerAMGSetCoarsenType(precond, 6);
      HYPRE_BoomerAMGSetOldDefault(precond);
      HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
      HYPRE_BoomerAMGSetNumSweeps(precond, 1);
      HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
      HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */
      /* Set the FlexGMRES preconditioner */
      HYPRE_FlexGMRESSetPrecond(linear_solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                          (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
      /* Now setup and solve! */
      HYPRE_ParCSRFlexGMRESSetup(linear_solver, parcsr_A, par_x, par_x);
   }
   else {
      printf ( "Use default linear solver.\n" );
   }

   //创建特征值求解器: 0 AMG, 1 PCG, 2 GMRES
   GCGE_SOLVER *hypre_solver = GCGE_HYPRE_Solver_Init_LinearSolverGivenByUser(parcsr_A, parcsr_B, linear_solver, solver_id, nev, argc, argv);   
   //一些参数的设置
   hypre_solver->para->ev_tol = 1e-8;
   hypre_solver->para->dirichlet_boundary = 0;
   hypre_solver->para->cg_max_it = 20;
   hypre_solver->para->ev_max_it = 150;
   hypre_solver->para->cg_type  = 2;

   //求解特征值问题
   GCGE_SOLVER_Solve(hypre_solver);  
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
   else if (solver_id == 2)
   {
      HYPRE_ParCSRFlexGMRESDestroy(linear_solver);
      HYPRE_BoomerAMGDestroy(precond);
   }
   else {

   }

   hypre_TFree(exact_eigenvalues);
   //释放矩阵空间
   HYPRE_IJMatrixDestroy(A);
   HYPRE_IJMatrixDestroy(B);
   HYPRE_IJVectorDestroy(x);

   /* Finalize MPI*/
   MPI_Finalize();
   return 0;
}
