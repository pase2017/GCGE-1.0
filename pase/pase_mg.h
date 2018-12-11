#ifndef _pase_mg_h_
#define _pase_mg_h_

#include "pase_ops.h"


#define pase_MultiGridDataAArray(data)  ((data)->A_array)
#define pase_MultiGridDataBArray(data)  ((data)->B_array)
#define pase_MultiGridDataPArray(data)  ((data)->P_array)
#define pase_MultiGridDataQArray(data)  ((data)->Q_array)
#define pase_MultiGridDataUArray(data)  ((data)->U_array)
#define pase_MultiGridDataFArray(data)  ((data)->F_array)

#define pase_MultiGridDataNumLevels(data) ((data)->num_levels)
      
typedef struct pase_MultiGrid_struct 
{
   PASE_Int           num_levels;
   void               **A_array;
   void               **B_array;
   void               **P_array;
   /* P0P1P2  P1P2  P2 */
   void               **Q_array;
   /* rhs and x */
   void               **U_array;
   void               **F_array;

   GCGE_OPS           *gcge_ops;
   
} pase_MultiGrid;
typedef struct pase_MultiGrid_struct *PASE_MultiGrid;

PASE_Int PASE_MultiGridCreate(PASE_MultiGrid* multi_grid, PASE_Int max_levels, 
   void *A, void *B, GCGE_OPS *pase_ops);

PASE_Int PASE_MultiGridDestroy(PASE_MultiGrid* multi_grid);

PASE_Int PASE_MultiGridFromItoJ(PASE_MultiGrid multi_grid, PASE_Int level_i, PASE_Int level_j, 
      PASE_Int num, void **pvx_i, void** pvx_j);


#endif
