
/*
 * PASE_Solver_Init:
 * 
 *                                  PASE_Convert
 *                                       +
 *                                   PASE_AMG
 *                                       +
 * PASE_Solver:1 APP--->GCGE_OPS--->PASE_MultiGrid--->PASE_MatVec
 *                                                        +
 *                                                   GCGE_OPS--->PASE_OPS--->GCGE_OPS
 *                                                                   |
 *                                                          GCGE_Solver_PASE_Init
 *
 *
 *             2 APP--->GCGE_OPS
 *                         |
 *                 GCGE_Solver_APP_Init
 *
 *             3 PASE_LinearSolver<---PASE_MultiGrid   
 * */
