/*
 * =====================================================================================
 *
 *       Filename:  gcge_orthogonal.c
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
#include <math.h>


#include "gcge_orthogonal.h"

//对V的所有列向量做关于矩阵B的正交化，如果B=NULL，那么进行L2正交化
//全部正交化，则start=0
//V是要正交化的向量组(GCGE_Vec),B是矩阵(GCGE_MATRIX)

void GCGE_Orthogonal(void **V, GCGE_INT start, GCGE_INT *end, void *B, GCGE_OPS *ops, 
      GCGE_ORTH_PARA *orth_para, GCGE_WORKSPACE *workspace)
{
   /* TODO */
#if 0
    GCGE_DOUBLE    t1 = GCGE_GetTime();
#endif
    //GCGE_Printf(MPI_COMM_WORLD, "line 458, start: %d, end: %d\n", start, *end);
    if(ops->Orthogonal)
    {
        //GCGE_Printf(MPI_COMM_WORLD, "line 473, start: %d, end: %d\n", start, *end);
        ops->Orthogonal(V, start, end, B, orth_para->orth_zero_tol, workspace->V_tmp);
    }
    else
    {
        GCGE_INT    i, j, n_nonzero = 0, n_zero = 0, reorth_time,
                    *Ind            = workspace->orth_ind,
                    print_orthzero  = orth_para->print_orth_zero,
                    max_reorth_time = orth_para->max_reorth_time;
        GCGE_DOUBLE vin, vout, tmp, dd, 
                    orth_zero_tol   = orth_para->orth_zero_tol,
                    reorth_tol      = orth_para->reorth_tol,
                    *inner_product  = workspace->subspace_dtmp;
        void        **V_tmp         = workspace->V_tmp;
        //L2正交化
        if(B == NULL)
        {
            //从V中start位置开始进行正交化
            for( i=start; i<(*end); i++ )
            {
                //如果i==0,说明start=0,对第一个向量直接进行单位化
                if(i == 0)
                {
                    dd = GCGE_VecNorm(V[0], ops);
                    if(dd > 10*orth_zero_tol)
                    {
                        //做归一化
                        ops->VecAxpby(0.0, V[0], 1.0/dd, V[0]);
                        Ind[0] = 0;
                        n_nonzero += 1;
                    }
                }
                else
                {
                    //vout = GCGE_VecNorm(V[i], ops);
                    reorth_time = 0;
                    //进行最多max_reorth_time次重正交化
                    do{
                        reorth_time += 1;
                        //先计算V[i]与前面的V[j]的局部内积，存放在inner_product中
                        for( j=0; j<start; j++ )
                        {
                            ops->VecLocalInnerProd(V[i], V[j], inner_product+j);
                        }
                        for( j=0; j<n_nonzero; j++ )
                        {
                            ops->VecLocalInnerProd(V[i], V[Ind[j]], inner_product+start+j);
                        }
                        ops->VecLocalInnerProd(V[i], V[i], inner_product+start+n_nonzero);
                        /* TODO */
//#if USE_MPI
#if 0
                        //共有start+n_nonzero+1个局部内积要发送
                        MPI_Allreduce(MPI_IN_PLACE, inner_product, start+n_nonzero+1, 
                                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
                        //现在inner_product中为V[i]与相应的V[j]的内积，下面减去相应的分量
                        for( j=0; j<start; j++ )
                        {
                            ops->VecAxpby(-inner_product[j], V[j], 1.0, V[i]);
                        }
                        for( j=0; j<n_nonzero; j++ )
                        {
                            ops->VecAxpby(-inner_product[start+j], V[Ind[j]], 1.0, V[i]);
                        }
                        //下面对V[i]进行归一化，首先计算V[i]的范数
                        vout = inner_product[start+n_nonzero];
                        vin = sqrt(vout);
                        for( j=0; j<start+n_nonzero; j++ )
                        {
                            vout -= inner_product[j]*inner_product[j];
                        }
                        vout = sqrt(vout);
                    }while((vout/vin < reorth_tol)&&(reorth_time < max_reorth_time));
                    //如果向量非0，进行单位化
                    if(vout > 10*orth_zero_tol)
                    {
                        ops->VecAxpby(0.0, V[i], 1.0/vout, V[i]);
                        Ind[n_nonzero++] = i;
                    }
                    else
                    {
                        if(print_orthzero == 1)
                        {
                            printf("In Orthogonal, there is a zero vector!, "
                                    "i = %d, start = %d, end = %d\n", i, start, *end);
                        }
                    }
                }
            }
        }
        else
        {
            //从V中start位置开始进行正交化
            //B正交化，先计算前start个向量的矩阵乘向量放到V_tmp中
            for( i=0; i<start; i++ )
            {
                ops->MatDotVec(B, V[i], V_tmp[i]);
            }
            for( i=start; i<(*end); i++ )
            {
                //如果i==0,说明start=0,对第一个向量直接进行单位化
                if(i == 0)
                {
                    dd = sqrt(GCGE_VecMatrixVec(V[0], B, V[0], V_tmp[0], ops));
                    if(dd > 10*orth_zero_tol)
                    {
                        ops->VecAxpby(0.0, V[0], 1.0/dd, V[0]);
                        Ind[n_nonzero++] = 0;
                        ops->MatDotVec(B, V[0], V_tmp[0]);
                    }
                }
                else
                {
                    //vout = sqrt(GCGE_VecMatrixVec(V[i], B, V[i], V_tmp[start+n_nonzero], ops));
                    reorth_time = 0;
                    //进行最多max_reorth_time次重正交化
                    do{
                        ops->MatDotVec(B, V[i], V_tmp[start+n_nonzero]);
                        reorth_time += 1;
                        //先计算V[i]与前面的V[j]的局部内积，存放在inner_product中
                        for( j=0; j<start+n_nonzero+1; j++ )
                        {
                            ops->VecLocalInnerProd(V[i], V_tmp[j], inner_product+j);
                        }
//#if USE_MPI
#if 0
                        //共有start+n_nonzero+1个局部内积要发送
                        MPI_Allreduce(MPI_IN_PLACE, inner_product, start+n_nonzero+1, 
                                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
                        //现在inner_product中为V[i]与相应的V[j]的内积，下面减去相应的分量
                        for( j=0; j<start; j++ )
                        {
                            ops->VecAxpby(-inner_product[j], V[j], 1.0, V[i]);
                        }
                        for( j=0; j<n_nonzero; j++ )
                        {
                            ops->VecAxpby(-inner_product[start+j], V[Ind[j]], 1.0, V[i]);
                        }
                        //下面对V[i]进行归一化，首先计算V[i]的范数
                        vout = inner_product[start+n_nonzero];
                        vin = sqrt(vout);
                        for( j=0; j<start+n_nonzero; j++ )
                        {
                            vout -= inner_product[j]*inner_product[j];
                        }
                        vout = sqrt(vout);
                    }while((vout/vin < reorth_tol)&&(reorth_time < max_reorth_time));
                    //如果向量非0，进行单位化
                    if(vout > 10*orth_zero_tol)
                    {
                        ops->VecAxpby(0.0, V[i], 1.0/vout, V[i]);
                        ops->MatDotVec(B, V[i], V_tmp[start+n_nonzero]);
                        Ind[n_nonzero++] = i;
                    }
                    else
                    {
                        if(print_orthzero == 1)
                        {
                            printf("In Orthogonal, there is a zero vector!, "
                                    "i = %d, start = %d, end = %d\n", i, start, *end);
                        }
                        //Nonzero_Vec[n_zero++] = V[i];
                    }
                }
            }
        }
        //接下来要把V的所有非零列向量存储在地址表格中靠前位置
        if(n_nonzero < (*end-start))
        {
            *end = start+n_nonzero;
            for( i=0; i<n_nonzero; i++ )
            {
                if(start+i != Ind[i])
                {
                    //V[start+i] = V[Ind[i]];
                    ops->VecAxpby(1.0, V[Ind[i]], 0.0, V[start+i]);
                }
            }
        }
    }
    //TODO 统计正交化时间
//    GCGE_DOUBLE t2 = GCGE_GetTime();
//    para->stat_para->PartTimeTotal->Orth_Time += t2-t1;
//    para->stat_para->PartTimeInOneIte->Orth_Time = t2-t1;
}

//进行部分的正交化, 对V中start位置中后的向量与前start的向量做正交化，
//同时V的start之后的向量自己也做正交化,ldV表示向量长度
void GCGE_OrthogonalSubspace(GCGE_DOUBLE *V, GCGE_INT start, GCGE_INT *end, GCGE_INT ldV, 
      GCGE_DOUBLE *B, GCGE_INT ldB, GCGE_ORTH_PARA *orth_para, GCGE_WORKSPACE *workspace)
{
    GCGE_INT    i, j, n_nonzero = 0, reorth_time = 0,
                print_orthzero  = orth_para->print_orth_zero,
                max_reorth_time = orth_para->max_reorth_time,
                *Ind            = workspace->orth_ind;
    GCGE_DOUBLE vin, vout, tmp, dd,
                orth_zero_tol = orth_para->orth_zero_tol, 
                reorth_tol    = orth_para->reorth_tol;

    if(B == NULL)
    {
        for( i=start; i<(*end); i++ )
        {
            if(i == 0)
            {
                dd = GCGE_VecNormSubspace(V, ldV);
                if(dd > 10*orth_zero_tol)
                {
                    GCGE_VecScaleSubspace(1.0/dd, V, ldV);
                    Ind[0] = 0;
                    n_nonzero = 1;
                }
            }
            else
            {
         
                vout = GCGE_VecNormSubspace(V+i*ldV, ldV);
                reorth_time = 0; 
                do{
                    vin = vout;
                    for(j = 0; j < start; j++)
                    {
                        tmp = GCGE_VecDotVecSubspace(V+j*ldV, V+i*ldV, ldV);
                        GCGE_VecAXPBYSubspace(-tmp, V+j*ldV, 1.0, V+i*ldV, ldV);
                    }
                    for(j = 0; j < n_nonzero; j++)
                    {
                        tmp = GCGE_VecDotVecSubspace(V+Ind[j]*ldV, V+i*ldV, ldV);
                        GCGE_VecAXPBYSubspace(-tmp, V+Ind[j]*ldV, 1.0, V+i*ldV, ldV);
                    }
                    vout = GCGE_VecNormSubspace(V+i*ldV, ldV);
                    reorth_time += 1;
                }while((vout/vin < reorth_tol)&&(reorth_time < max_reorth_time));
         
                if(vout > 10*orth_zero_tol)
                {
                    GCGE_VecScaleSubspace(1.0/vout, V+i*ldV, ldV);
                    Ind[n_nonzero++] = i;
                }
                else
                {
                    if(print_orthzero == 1)
                    {
                        printf("In OrthogonalSubspace, there is a zero vector!"
                               "i = %d, start = %d, end: %d\n", i, start, *end);
                    }
                }
            }

        }
    }

    if(n_nonzero < (*end-start))
    {
        *end = start+n_nonzero;
        for( i=0; i<n_nonzero; i++ )
        {
            if((start+i) != Ind[i])
            {
                //memcpy(V+(start+i)*ldV, V+Ind[i]*ldV, ldV*sizeof(GCGE_DOUBLE));
                GCGE_VecCopySubspace(V+Ind[i]*ldV, V+(start+i)*ldV, ldV);
            }
        }
    }
}
