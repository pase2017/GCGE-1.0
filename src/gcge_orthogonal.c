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
 * 2018-10-03 by Hehu Xie
 *  改了一下正交化的选择：精度低的时候选择block的形式，当精度高了之后选择单个向量的形式,
 * 同时GCGE_Orthogonal函数的para-Orth_para 改成了para,为了得到sum_res的值.
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "gcge_orthogonal.h"


/**
 * @brief 
 *
 * 正交化V, 从start到end, 且返回end, 因为正交化过程中会去掉0向量
 *
 * V_tmp只需要一个，记录当前要正交的向量V[current]与B的乘积
 * d_tmp记录当前要正交的向量V[current]与之前向量以及与自己的B内积, 长度为输入的end
 *
 * 过程中需要做重正交化，当向量的模改变太大(此时，会包含太多以前的分量)
 *
 * 若向量的模太小，则与*end向量交换指针，并(*end)--
 *
 * 最后做一次归一化
 *
 * 注意：大写表示向量组，小写表示单向量
 *
 * for current = start ; current < (*end); ++current
 *     计算V_tmp
 *     reorth_count = 0;
 *     do{
 *     MatDotVec(B, V[current], V_tmp)
 *     计算d_tmp
 *
 *     mv_s = {0, 0}
 *     mv_e = {current+1, 1}
 *     MultiVecInnerProd(V, &V_tmp, d_tmp, "ns", d_tmp, mv_s, mv_e, current+1, ops )
 *
 *    // V[current] = V[current] - sum_{i=0}^{current-1} d_tmp[i] V[i]
 *    // d_tmp[current] = d_tmp[current] - sum_{i=0}^{current-1} d_tmp[i]^2
 *    // V[current] = V[current] / d_tmp[current]
 * 
 *     if(current != 0)
 *     {
 *       mv_s = {0, 0}
 *       mv_e = {current, 1}
 *       MultiVecLinearComb (V, &V_tmp, mv_s, mv_e, d_tmp, current, ops)
 *       VecAxpby(-1, V_tmp, 1, V[current]);
 *     }
 *
 *     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *     记 vi0 为正交化前的 v_i, vi1 为去掉前面分量后的 v_i
 *     那么有 vi1 = vi0 - \sum_{j=0}^{i-1} (vi0, v_j) * v_j
 *     要计算vi1的范数，我们有
 *     (vi1, vi1) = (vi0 - \sum_{j=0}^{i-1} (vi0, v_j) * v_j, vi0 - \sum_{j=0}^{i-1} (vi0, v_j) * v_j)
 *                = (vi0, vi0) - \sum_{j=0}^{i-1} (vi0, v_j)^2
 *     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *
 *     vin  = sqrt(d_tmp[current])
 *     d_tmp[current] = d_tmp[current] - sum_{i=0}^{current-1} d_tmp[i]^2
 *     vout = sqrt(d_tmp[current])
 *
 *     ratio = d_tmp[current] / vin
 *
 *     reorth_count++
 *     
 *     }while(ratio < orth_para->reorth_tol && reorth_count < orth_para->max_reorth_count)
 *
 *     if (d_tmp[current] > orth_para->orth_zero_tol)
 *     {
 *         VecAxpby(0, V_tmp, 1/d_tmp[current], V[current])
 *     }
 *     else 
 *     {
 *         mv_s = {current, *end-1}
 *         mv_e = {current+1, *end}
 *         MultiVecSwap(V, V, mv_s, mv_e, ops);
 *         (*end)--;
 *         current--;
 *     }
 * end
 *
 * @param V
 * @param start 从0 到 start-1已经正交归一
 * @param end   输入要正交的向量个数， 返回被正交化的向量个数
 * @param B
 * @param ops
 * @param orth_para
 * @param V_tmp 长度1
 * @param d_tmp 长度×end
 */
 //对V的某些列进行正交化，start: 表示开始的位置，end[0]表示之前已经做好正交化向量的最后位置，
 //end[1]表示需要正交化向量的最终位置，矩阵B表示进行内积计算的矩阵，orth_para: 定义正交化方式的参数
 //V_tmp: 用来存储正交化过程中的临时向量, d_tmp: 用来存储当前向量与之前已经正交化向量的内积 
void GCGE_Orthogonal(void **V, GCGE_INT start, GCGE_INT *end, 
      void *B, GCGE_OPS *ops, GCGE_PARA *para, void **V_tmp, GCGE_DOUBLE *d_tmp)
{
    GCGE_ORTH_PARA *orth_para = para->orth_para; 
    GCGE_INT    current = 0; //表示for循环中当前正在进行正交化的向量编号
    GCGE_INT    reorth_count = 0; //统计(重)正交化次数
    //mv_s: 用来记录正交化的初始和终止向量的位置
    //mv_start
    GCGE_INT    mv_s[2] = { 0, 0 };
    //mv_end: 
    GCGE_INT    mv_e[2] = { 0, 0 };

    GCGE_DOUBLE vin = 0.0, vout = 0.0, norm_value = 0.0, ipv = 0.0; //ipv: inner product value
    GCGE_DOUBLE ratio = 1.0, sum_res = para->sum_res, 
                criterion_tol = orth_para->criterion_tol;
    void        *vec_current; //当前正在进行正交化的向量
    void        *vec_tmp, *vec_tmp2;
    GCGE_INT    j = 0;
    //ops->GetVecFromMultiVec(V_tmp, 1, &vec_tmp2);
    //end[0]: 表示需要正交化向量的最终位置
    for(current = start; current < (*end); ++current)
    {
        reorth_count = 0;

        //下面是严格的单个向量重正交化的方式(数值稳定性比较好)
        //这里根据目前的残差的量级来进行自动选择
        //printf("sum_res= %e\n",sum_res);
        if((sum_res >0.0)&&(sum_res<criterion_tol))
        {
            //printf("do the column version!\n");
            ops->GetVecFromMultiVec(V, current, &vec_current); //vec_current = V[current]
            if(B == NULL)
            {
                ops->VecInnerProd(vec_current, vec_current, &norm_value);      
                vout = sqrt(norm_value);    
            }
            else
            {
                //vec_tmp2 = B*vec_current
                ops->GetVecFromMultiVec(V_tmp, 1, &vec_tmp2); //vec_tmp2 = V_tmp[1]
                ops->MatDotVec(B, vec_current, vec_tmp2); 
                ops->VecInnerProd(vec_current, vec_tmp2, &norm_value);
                vout = sqrt(norm_value);
                ops->RestoreVecForMultiVec(V_tmp, 1, &vec_tmp2);
            }
            do{
                // B != NULL, 计算 vec_tmp = B * vec_current
                // B == NULL, 计算 vec_tmp = vec_current
                //ops->GetVecFromMultiVec(V_tmp, 0, &vec_tmp); //vec_tmp = V_tmp[0] 
                //ops->GetVecFromMultiVec(V_tmp, 1, &vec_tmp2); //vec_tmp2 = V_tmp[1]
                for (j=0; j<current; j++)
                {
                    ops->GetVecFromMultiVec(V, j, &vec_tmp); //vectmp = V(:,j)
                    //下面做内积
                    if(B == NULL)
                    {
                        ops->VecInnerProd(vec_current, vec_tmp, &ipv);
                    }
                    else
                    {
                        //vec_tmp2 = B*vec_current
                        ops->GetVecFromMultiVec(V_tmp, 1, &vec_tmp2); //vec_tmp2 = V_tmp[1]
                        ops->MatDotVec(B, vec_current, vec_tmp2); 
                        ops->VecInnerProd(vec_tmp, vec_tmp2, &ipv);        
                        ops->RestoreVecForMultiVec(V_tmp, 1, &vec_tmp2);        
                    }//end for if
                    //接下来做： V(:,current) = V(:,current) - ipv*V(:,j)
                    ops->VecAxpby(-ipv, vec_tmp, 1.0, vec_current);
                    //ops->GetVecFromMultiVec(V, j, &vec_tmp); //vectmp = V(:,j)
                    ops->RestoreVecForMultiVec(V, j, &vec_tmp);
                }//end for j
                vin = vout;       
                if(B == NULL)
                {
                    ops->VecInnerProd(vec_current, vec_current, &norm_value);
                    //vout = sqrt(norm_value);              
                }
                else
                {
                    //vec_tmp2 = B*vec_current
                    ops->GetVecFromMultiVec(V_tmp, 1, &vec_tmp2); //vec_tmp2 = V_tmp[1]
                    ops->MatDotVec(B, vec_current, vec_tmp2); 
                    ops->VecInnerProd(vec_current, vec_tmp2, &norm_value);
                    //vout = sqrt(norm_value);
                    ops->RestoreVecForMultiVec(V_tmp, 1, &vec_tmp2);          
                }
                //ops->VecInnerProd(vec_current, vec_current, &norm_value);
                vout = sqrt(norm_value);

                ratio = vout/vin;
                reorth_count;     
            }while((ratio < orth_para->reorth_tol) && (reorth_count < orth_para->max_reorth_time) 
                    && (vout > orth_para->orth_zero_tol) );
            //ops->RestoreVecForMultiVec(V_tmp, 0, &vec_tmp);
            //ops->RestoreVecForMultiVec(V_tmp, 1, &vec_tmp2);
        }
        else
        {
            printf("do the block version!\n");
            do{
                // B != NULL, 计算 vec_tmp = B * vec_current
                // B == NULL, 计算 vec_tmp = vec_current
                ops->GetVecFromMultiVec(V_tmp, 0, &vec_tmp);
                ops->GetVecFromMultiVec(V, current, &vec_current);
                if(B == NULL)
                {
                    ops->VecAxpby(1.0, vec_current, 0.0, vec_tmp);
                }
                else
                {
                    ops->MatDotVec(B, vec_current, vec_tmp);
                }
                ops->RestoreVecForMultiVec(V_tmp, 0, &vec_tmp);
                ops->RestoreVecForMultiVec(V, current, &vec_current);

                // 计算 d_tmp = V^T * vec_tmp
                // 即为 d_tmp[j] = (v_j, vec_current)_B
                mv_s[0] = 0;
                mv_e[0] = current+1;
                mv_s[1] = 0;
                mv_e[1] = 1;
                ops->MultiVecInnerProd(V, V_tmp, d_tmp, "ns", mv_s, mv_e, current+1, ops);

                // 如果当前计算的不是向量组中的第一个向量，那么需要减去前面的分量
                if(current != 0)
                {
                    // 计算线性组合 V_tmp[0] = V * d_tmp
                    // 即为 V_tmp[0] = \sum_{j=0}^{current-1} d_tmp[j] * v_j
                    mv_s[0] = 0;
                    mv_e[0] = current;
                    mv_s[1] = 0;
                    mv_e[1] = 1;
                    ops->MultiVecLinearComb(V, V_tmp, mv_s, mv_e, d_tmp, current, NULL, -1, ops);

                    // 计算 vec_current = vec_current - V_tmp[0]
                    ops->GetVecFromMultiVec(V_tmp, 0, &vec_tmp);
                    ops->GetVecFromMultiVec(V, current, &vec_current);
                    ops->VecAxpby(-1.0, vec_tmp, 1.0, vec_current);
                    ops->RestoreVecForMultiVec(V_tmp, 0, &vec_tmp);
                    ops->RestoreVecForMultiVec(V, current, &vec_current);
                }

                // d_tmp[current] 即为 (v_current_in, vec_current_in)_B
                // vin为vec_current去掉前面分量前的B范数
                vin = sqrt(d_tmp[current]);
                /* 记 vi0 为正交化前的 v_i, vi1 为去掉前面分量后的 v_i
                 *  那么有 vi1 = vi0 - \sum_{j=0}^{i-1} (vi0, v_j) * v_j
                 *  要计算vi1的范数，我们有
                 *  (vi1, vi1) = (vi0 - \sum_{j=0}^{i-1} (vi0, v_j) * v_j, 
                 *                vi0 - \sum_{j=0}^{i-1} (vi0, v_j) * v_j)
                 *             = (vi0, vi0) - \sum_{j=0}^{i-1} (vi0, v_j)^2 */
                d_tmp[current] = d_tmp[current] - GCGE_VecDotVecSubspace(d_tmp, d_tmp, current);
                // ratio = vout / vin
                vout = sqrt(d_tmp[current]);
                ratio = vout/vin;
                reorth_count++;
                //如果vec_current范数几乎没有变小，或者已经达到了最大重正交化次数
                //或者vec_current范数已经比判定向量为0的阈值还小
                //就不做重正交化

            }while((ratio < orth_para->reorth_tol) && (reorth_count < orth_para->max_reorth_time)
                    && (vout > orth_para->orth_zero_tol));
        }//end if (sum_res)


        //下面看看目前处理的向量是否是一个零向量    

        if(vout > orth_para->orth_zero_tol)
        {
            //如果vec_current不为0，就做归一化
            //ops->GetVecFromMultiVec(V, current, &vec_current);
            ops->VecAxpby(0.0, vec_current, 1/vout, vec_current);     
            ops->RestoreVecForMultiVec(V, current, &vec_current);
            //ops->RestoreVecForMultiVec(V_tmp, 0, &vec_tmp);
            //ops->RestoreVecForMultiVec(V_tmp, 1, &vec_tmp2);
        }
        else 
        {
            //如果vec_current为0，就把此向量与最后一个向量进行指针交换
            //同时(*end)--,把总向量的个数减1
            //同时current--,把当前向量的编号减1,循环到下次仍计算当前编号的向量(即原end)
            mv_s[0] = current;
            mv_s[1] = *end-1;
            mv_e[0] = current+1;
            mv_e[1] = *end;
            ops->MultiVecSwap(V, V, mv_s, mv_e, ops);
            (*end)--;
            current--;
            if(orth_para->print_orth_zero == 1)
            {
                printf("In Orthogonal, there is a zero vector!, "
                        "current = %d, start = %d, end = %d\n", current, start, *end);
            }//end if(orth_para->print_orth_zero == 1)      
        }//end if(vout > orth_para->orth_zero_tol)    

    }//end for(current = start; current < (*end); ++current)
}//end of this subprogram

/**
 * @brief 
 *
 * 正交化V, 从start到end, 且返回end, 因为正交化过程中会去掉0向量
 *
 * V_tmp只需要一个，记录当前要正交的向量V[current]与B的乘积
 * d_tmp记录当前要正交的向量V[current]与之前向量以及与自己的B内积, 长度为输入的end
 *
 * 过程中需要做重正交化，当向量的模改变太大(此时，会包含太多以前的分量)
 *
 * 若向量的模太小，则与*end向量交换指针，并(*end)--
 *
 * 最后做一次归一化
 *
 *
 * for current = start ; current < (*end); ++current
 *     计算V_tmp
 *     reorth_count = 0;
 *     do{
 *
 *    // V_tmp = B * V[current]
 *    // d_tmp = V^T * V_tmp
 *    // V[current] = V[current] - sum_{i=0}^{current-1} d_tmp[i] V[i]
 *    // d_tmp[current] = d_tmp[current] - sum_{i=0}^{current-1} d_tmp[i]^2
 *    // V[current] = V[current] / d_tmp[current]
 *
 *     ratio = d_tmp[current] / vin
 *
 *     reorth_count++
 *     
       }while(ratio < orth_para->reorth_tol && reorth_count < orth_para->max_reorth_time && vout > orth_para->orth_zero_tol);
 *
 *     if (d_tmp[current] > orth_para->orth_zero_tol)
 *         V[current] = V[current]/vout
 *     else 
 *         V[current] = V[*end-1]
 *         (*end)--;
 *         current--;
 *     end
 * end
 *
 * @param V
 * @param start 从0 到 start-1已经正交归一
 * @param end   输入要正交的向量个数， 返回被正交化的向量个数
 * @param B
 * @param ops
 * @param orth_para
 * @param V_tmp 长度1
 * @param d_tmp 长度×end
 */
void GCGE_OrthogonalSubspace(double *V, GCGE_INT ldV, GCGE_INT start, GCGE_INT *end, 
      void *B, GCGE_INT ldB, GCGE_ORTH_PARA *orth_para)
{
    GCGE_INT    current = 0; //当前进行正交化操作的向量编号
    GCGE_INT    reorth_count = 0;
    GCGE_DOUBLE vin = 0.0;
    GCGE_DOUBLE vout = 0.0;
    GCGE_DOUBLE ratio = 1.0;
    GCGE_INT    idx = 0;
    GCGE_DOUBLE ip = 1.0; //inner product
    for(current = start; current < (*end); ++current)
    {
        // 计算当前向量的范数
        vout = GCGE_VecNormSubspace(V+current*ldV, ldV);
        if(current != 0)
        {
            reorth_count = 0;
            do{
                vin = vout;
                for(idx = 0; idx < current; idx++)
                {
                    //修正的Gram-Schmidt正交化方法
                    //计算ip = (v_idx, v_current)
                    ip = GCGE_VecDotVecSubspace(V+idx*ldV, V+current*ldV, ldV);
                    //计算v_current = v_current - ip * v_idx
                    GCGE_VecAXPBYSubspace(-ip, V+idx*ldV, 1.0, V+current*ldV, ldV);
                }
                // 计算当前向量的范数
                vout = GCGE_VecNormSubspace(V+current*ldV, ldV);

                ratio = vout/vin;

                reorth_count++;

            }while(ratio < orth_para->reorth_tol && reorth_count < orth_para->max_reorth_time && vout > orth_para->orth_zero_tol);
        }

        if(vout > orth_para->orth_zero_tol)
        {
            //如果v_current不为0，就做归一化
            GCGE_VecScaleSubspace(1.0/vout, V+current*ldV, ldV);
        }//做完归一化的情况，下面考虑是 0 向量的情况
        else 
        {
            //如果v_current为0，就把最后一个向量拷贝到当前向量的位置
            //如果本向量就是最后一个向量，那就不拷贝(考虑如果是最后一个直接跳出)
            //同时(*end)--,把总向量的个数减1
            //同时current--,把当前向量的编号减1,循环到下次仍计算当前编号的向量(即原end)
            if(current != *end-1)
            {
                   //下面这个函数的方式应该是copy（from ,to， OK！
                GCGE_VecCopySubspace(V+(*end -1)*ldV, V+current*ldV, ldV);
            }//end if (current != *end-1)
            (*end)--;
            current--;
            //如果这个时候需要输出提醒信息，就输出
            if(orth_para->print_orth_zero == 1)
            {
                printf("In OrthogonalSubspace, there is a zero vector!, "
                        "current = %d, start = %d, end = %d\n", current, start, *end);
            }//end if (orth_para->print_orth_zero == 1)
        }//end if (vout > orth_para->orth_zero_tol)
    }//end for current to the end
}//end for this subprogram

