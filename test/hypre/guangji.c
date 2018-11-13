#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "math.h"


void Read_dmatcsr(const char *filename, 
	          int*  nr, int*  nc, int*     nn,
		  int** ia, int** ja, double** va)
{
    FILE *file = fopen(filename, "r");
    if(!file)
    {
        printf("\nError: Cannot open %s!\n", filename);
	exit(-1);
    }

    int i = 0;
    while(i<1 && EOF!=fscanf(file, "%d\n",   nr)) i++;
    i = 0;
    while(i<1 && EOF!=fscanf(file, "%d\n",   nc)) i++;
    i = 0;
    while(i<1 && EOF!=fscanf(file, "%d\n\n", nn)) i++;
    i = 0;

    *ia = (int*)   malloc((*nr+1)*sizeof(int));
    *ja = (int*)   malloc((*nn)  *sizeof(int));
    *va = (double*)malloc((*nn)  *sizeof(double));

    while(i<*nr+1 && EOF!=fscanf(file, "%d\n", *ia+i))  i++;
    i = 0;
    while(i<1     && EOF!=fscanf(file, "\n"))           i++;
    i = 0;
    while(i<*nn   && EOF!=fscanf(file, "%d\n", *ja+i))  i++;
    i = 0;
    while(i<1     && EOF!=fscanf(file, "\n"))           i++;
    i = 0;
    while(i<*nn   && EOF!=fscanf(file, "%lf\n", *va+i)) i++;

    fclose(file);
}

HYPRE_Int CreateGuangJiEigenProblem(MPI_Comm comm, 
   HYPRE_IJMatrix* A, HYPRE_IJMatrix* B, HYPRE_IJMatrix* Asub, 
   HYPRE_IJVector* x, HYPRE_IJVector* b )
{
   int i, N;
   HYPRE_Int myid, num_procs;
   HYPRE_Int local_size, extra, ilower, iupper; 

   MPI_Comm_rank(comm, &myid);
   MPI_Comm_size(comm, &num_procs);

   HYPRE_Int    A_nr,     A_nc,     A_nn;
   HYPRE_Int   *A_ia,    *A_ja;
   HYPRE_Real  *A_va;
   HYPRE_Int    B_nr,     B_nc,     B_nn;
   HYPRE_Int   *B_ia,    *B_ja;
   HYPRE_Real  *B_va;
   HYPRE_Int    Asub_nr,  Asub_nc,  Asub_nn;
   HYPRE_Int   *Asub_ia, *Asub_ja;
   HYPRE_Real  *Asub_va;

   if(myid == 0) {
     printf("Reading stiff matrixing....\n");
   }
   Read_dmatcsr("/share/home/hhxie/MATRIX/guangji_273w/guangji.stif.dmatcsr", &A_nr, &A_nc, &A_nn, &A_ia, &A_ja, &A_va);
   N = A_nr;

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

   HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, A);
   HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, B);
   HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, Asub);

   HYPRE_IJMatrixSetObjectType(*A, HYPRE_PARCSR);
   HYPRE_IJMatrixSetObjectType(*B, HYPRE_PARCSR);
   HYPRE_IJMatrixSetObjectType(*Asub, HYPRE_PARCSR);

   /* Initialize before setting coefficients */
   HYPRE_IJMatrixInitialize(*A);
   HYPRE_IJMatrixInitialize(*B);
   HYPRE_IJMatrixInitialize(*Asub);

   if(myid == 0) {
     printf("Set values for stiff matrix......\n");
   }
   {
      int nnz;
      for (i = ilower; i <= iupper; i++)
      {
         nnz = A_ia[i+1]-A_ia[i];
         /* Set the values for row i */
         HYPRE_IJMatrixSetValues(*A, 1, &nnz, &i, A_ja+A_ia[i], A_va+A_ia[i]);
      }
   }
   free(A_ia);  free(A_ja);  free(A_va);
   if(myid == 0) {
     printf("Reading mass matrixing....\n");
   }
   Read_dmatcsr("/share/home/hhxie/MATRIX/guangji_273w/guangji.mass.dmatcsr", &B_nr, &B_nc, &B_nn, &B_ia, &B_ja, &B_va);
   if(myid == 0) {
     printf("Set values for mass matrix......\n");
   }
   {
      int nnz;
      for (i = ilower; i <= iupper; i++)
      {
         nnz = B_ia[i+1]-B_ia[i];
         /* Set the values for row i */
         HYPRE_IJMatrixSetValues(*B, 1, &nnz, &i, B_ja+B_ia[i], B_va+B_ia[i]);
      }
   }
   free(B_ia);  free(B_ja);  free(B_va);
   if(myid == 0) {
     printf("Reading Asub matrixing....\n");
   }
   Read_dmatcsr("/share/home/hhxie/MATRIX/guangji_273w/Asubmatrix.dmatcsr", &Asub_nr, &Asub_nc, &Asub_nn, &Asub_ia, &Asub_ja, &Asub_va);
   if(myid == 0) {
     printf("Set values for Asub matrix......\n");
   }
   {
      int nnz;
      for (i = ilower; i <= iupper; i++)
      {
         nnz = Asub_ia[i+1]-Asub_ia[i];
         /* Set the values for row i */
         HYPRE_IJMatrixSetValues(*Asub, 1, &nnz, &i, Asub_ja+Asub_ia[i], Asub_va+Asub_ia[i]);
      }
   }
   free(Asub_ia);
   free(Asub_ja);
   free(Asub_va);
   if(myid == 0) {
     printf("done......\n");
   }
   printf("N = %d, iupper = %d, ilower = %d, local_size = %d\n", N, iupper, ilower, local_size);

   /* Assemble after setting the coefficients */
   HYPRE_IJMatrixAssemble(*A);
   HYPRE_IJMatrixAssemble(*B);
   HYPRE_IJMatrixAssemble(*Asub);

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
HYPRE_Int DestroyGuangJiEigenProblem(
   HYPRE_IJMatrix A, HYPRE_IJMatrix B, HYPRE_IJMatrix Asub, 
   HYPRE_IJVector x, HYPRE_IJVector b)
{	   
   HYPRE_IJMatrixDestroy(A);
   HYPRE_IJMatrixDestroy(B);
   HYPRE_IJMatrixDestroy(Asub);
   HYPRE_IJVectorDestroy(x);
   HYPRE_IJVectorDestroy(b);

   return 0;  
}

