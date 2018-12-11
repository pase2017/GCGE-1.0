
/*
 * GCGE将接口传给PASE, PASE生成自己的矩阵向量操作
 * PASE再以此调用GCGE, 即将其矩阵向量操作再赋给GCGE
 *
 * SLEPC->GCGE-------->PASE->GCGE
 * */
