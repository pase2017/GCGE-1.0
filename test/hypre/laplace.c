#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "math.h"

static int cmp( const void *a ,  const void *b )
{   return *(double *)a > *(double *)b ? 1 : -1; }

HYPRE_Int CreateLaplaceEigenProblem(MPI_Comm comm, HYPRE_Int dim, 
   HYPRE_Real** exact_eigenvalues, HYPRE_Int block_size, 
   HYPRE_Int n,
   HYPRE_IJMatrix* A, HYPRE_IJMatrix* B, 
   HYPRE_IJVector* x, HYPRE_IJVector* b )
{
   int i, k, N;
   HYPRE_Int myid, num_procs;
   HYPRE_Int local_size, extra, ilower, iupper; 
   HYPRE_Real h, h2;

   MPI_Comm_rank(comm, &myid);
   MPI_Comm_size(comm, &num_procs);

   /* for 2D */
   /* Preliminaries: want at least one processor per row */
   if (n*n < num_procs) n = sqrt(num_procs) + 1;
   N = n*n; /* global number of rows */
   h = 1.0/(n+1); /* mesh size*/
   h2 = h*h;

   /*----------------------- Laplace精确特征值 ---------------------*/
   {
      int tmp_nn = (int) sqrt(block_size) + 3;
      *exact_eigenvalues = hypre_CTAlloc (HYPRE_Real, tmp_nn*tmp_nn);
      for (i = 0; i < tmp_nn; ++i) 
      {
	 for (k = 0; k < tmp_nn; ++k) 
	 {
	    (*exact_eigenvalues)[i*tmp_nn+k] = ( 4*sin( (i+1)*M_PI/(2*(n+1)) )*sin( (i+1)*M_PI/(2*(n+1)) ) 
		  + 4*sin( (k+1)*M_PI/(2*(n+1)) )*sin( (k+1)*M_PI/(2*(n+1)) ) );
	 }
      }
      qsort((*exact_eigenvalues), tmp_nn*tmp_nn, sizeof(double), cmp);
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
   HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, A);
   HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, B);

   /* Choose a parallel csr format storage (see the User's Manual) */
   HYPRE_IJMatrixSetObjectType(*A, HYPRE_PARCSR);
   HYPRE_IJMatrixSetObjectType(*B, HYPRE_PARCSR);

   /* Initialize before setting coefficients */
   HYPRE_IJMatrixInitialize(*A);
   HYPRE_IJMatrixInitialize(*B);

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
	 HYPRE_IJMatrixSetValues(*A, 1, &nnz, &i, cols, values);
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
//	 values[0] = h2;
	 values[0] = 1.0;
	 /* Set the values for row i */
	 HYPRE_IJMatrixSetValues(*B, 1, &nnz, &i, cols, values);
      }
   }
   /* Assemble after setting the coefficients */
   HYPRE_IJMatrixAssemble(*A);
   HYPRE_IJMatrixAssemble(*B);

   /* Create sample rhs and solution vectors */
   HYPRE_IJVectorCreate(comm, ilower, iupper, x);
   HYPRE_IJVectorSetObjectType(*x, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(*x);
   HYPRE_IJVectorAssemble(*x);  

   HYPRE_IJVectorCreate(comm, ilower, iupper, b);
   HYPRE_IJVectorSetObjectType(*b, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(*b);
   HYPRE_IJVectorAssemble(*b); 

   return 0;
}   

HYPRE_Int DestroyLaplaceEigenProblem(HYPRE_Real* exact_eigenvalues, 
   HYPRE_IJMatrix A, HYPRE_IJMatrix B, 
   HYPRE_IJVector x, HYPRE_IJVector b)
{	   
   hypre_TFree(exact_eigenvalues);
   HYPRE_IJMatrixDestroy(A);
   HYPRE_IJMatrixDestroy(B);
   HYPRE_IJVectorDestroy(x);
   HYPRE_IJVectorDestroy(b);

   return 0;  
}
