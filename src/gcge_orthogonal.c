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
 * 改了一下正交化的选择：精度低的时候选择block的形式，当精度高了之后选择单个向量的形式,
 * 同时GCGE_Orthonormalization函数的para-Orth_para 改成了para,为了得到sum_res的值.
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
void GCGE_Orthonormalization(void **V, GCGE_INT start, GCGE_INT *end, 
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
        if((sum_res > 0.0)&&(sum_res < criterion_tol))
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
                }//end if(B == NULL)
                //ops->VecInnerProd(vec_current, vec_current, &norm_value);
                vout = sqrt(norm_value);
                ratio = vout/vin;
                reorth_count;      
            }while((ratio < orth_para->reorth_tol) && (reorth_count < orth_para->max_reorth_time) 
                    && (vout > orth_para->orth_zero_tol) );
            //ops->RestoreVecForMultiVec(V_tmp, 0, &vec_tmp);
            //ops->RestoreVecForMultiVec(V_tmp, 1, &vec_tmp2);       
            ops->RestoreVecForMultiVec(V, current, &vec_current); //vec_current = V[current]
        }//结束了单个向量的正交化方式
        else
        {
            //printf("do the block version!\n");
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
                }//end if(current != 0)        

                // d_tmp[current] 即为 (v_current_in, vec_current_in)_B
                // vin为vec_current去掉前面分量前的B范数
                vin = sqrt(d_tmp[current]);
                /* 记 vi0 为正交化前的 v_i, vi1 为去掉前面分量后的 v_i
                 *  那么有 vi1 = vi0 - \sum_{j=0}^{i-1} (vi0, v_j) * v_j
                 *  要计算vi1的范数，我们有
                 *  (vi1, vi1) = (vi0 - \sum_{j=0}^{i-1} (vi0, v_j) * v_j, 
                 *                vi0 - \sum_{j=0}^{i-1} (vi0, v_j) * v_j)
                 *             = (vi0, vi0) - \sum_{j=0}^{i-1} (vi0, v_j)^2 */
                d_tmp[current] = d_tmp[current] - GCGE_ArrayDotArrayInSubspace(d_tmp, d_tmp, current);
                // ratio = vout / vin
                vout = sqrt(d_tmp[current]);
                ratio = vout/vin;
                reorth_count ++;
                //如果vec_current范数几乎没有变小，或者已经达到了最大重正交化次数
                //或者vec_current范数已经比判定向量为0的阈值还小
                //就不做重正交化        
            }while((ratio < orth_para->reorth_tol) && (reorth_count < orth_para->max_reorth_time)
                    && (vout > orth_para->orth_zero_tol));      
        }//end if (sum_res)做完块的形式了

        //下面看看目前处理的向量是否是一个零向量    

        if(vout > orth_para->orth_zero_tol)
        {
            //如果vec_current不为0，就做归一化
            ops->GetVecFromMultiVec(V, current, &vec_current);
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
                GCGE_Printf("In Orthonormalization, there is a zero vector!, "
                        "current = %d, start = %d, end = %d\n", current, start, *end);    
            }//end if(orth_para->print_orth_zero == 1)            
        }//end if(vout > orth_para->orth_zero_tol)      
    }//end for(current = start; current < (*end); ++current)    
}//end of this subprogram


//对V的某些列进行正交化，start: 表示开始的位置，end[0]表示需要正交化向量的最终位置(也返回这个最终位置)，
//矩阵B表示进行内积计算的矩阵，orth_para: 定义正交化方式的参数
//evec, V_tmp: 用来存储正交化过程中 [B*X, B*P] 的临时向量, 
//d_tmp: 用来存储待正交化的向量与之前已经正交化的一个向量的内积 
/** 
 *  V:         input,        要正交化的向量组
 *  start:     input,        要正交化的向量组在V中的起始位置
 *  end:       input|output, 要正交化的向量组在V中的终止位置
 *  B:         input,        正交化矩阵, 可以是NULL
 *  ops:       input,        操作
 *  para:      input,        参数
 *  workspace: input,        工作空间
 *
 *  块正交化过程：
 *  
 *  V = [V1, V2, V3], V中有三个部分, 其中V1,V2部分已经正交, 
 *  V1中有 size_V1(这里是nev) 个向量, V2中有 size_V2(这里是dim_xp-nev) 个向量
 *  V3中有 size_V3(这里是w_length=end-start) 个向量
 *
 *  临时存储空间使用workspace中的 evec(nev个), V_tmp(dim_xp-nev个), CG_p(1个)
 *  1. 首先计算 evec  = B * V1
 *              V_tmp = B * V2
 *
 *  下面的过程进行两次(相当于做一次重正交化)
 *  2. 去掉 V3 中 V1,V2 方向的分量
 *     for current = start_V1:end_V2
 *         计算 d_tmp = (V3, V[current])
 *         for V3_current = start_V3:end_V3
 *             计算 V[V3_current] = V[V3_current] - d_tmp[V3_current-start_V3] * V[current]
 *         end
 *     end
 *
 *  3. V3 自身做正交化 
 *     for current = start_V3:end_V3
 *         计算 d_tmp = (V[current:end_V3], V[current])
 *         norm = sqrt(d_tmp[0])
 *         如果norm > orth_zero_tol
 *            for V3_current = current:end_V3
 *                计算 V[V3_current] = V[V3_current] - d_tmp[V3_current-i]/norm * V[current]
 *            end
 *         否则，判定向量V[V3_current]为0
 *     end
 *
 */
void GCGE_BOrthonormalization(void **V, GCGE_INT start, GCGE_INT *end, 
                      void *B, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
    GCGE_ORTH_PARA *orth_para = para->orth_para; 
    GCGE_INT        current = 0, w_length = 0; //表示for循环中当前正在进行正交化的向量编号
    GCGE_INT        reorth_count = 0; //统计(重)正交化次数
    //mv_s: 用来记录正交化的初始和终止向量的位置
    //mv_start
    GCGE_INT    mv_s[2] = { 0, 0 };
    //mv_end: 
    GCGE_INT    mv_e[2] = { 0, 0 };

    GCGE_INT dim_x = workspace->dim_x, dim_xp = workspace->dim_xp, nev = para->nev;

    GCGE_DOUBLE norm_value = 0.0;  //存储当前向量的范数
    GCGE_DOUBLE *d_tmp = workspace->subspace_dtmp;
    void        *vec_current; //当前正在进行正交化的向量
    void        *vec_tmp;
    //CG_p[0] = B*W(:,j): 将用来存储临时向量， evec：用来存储  B*V（：，0：nev-1）， V_tmp = B*V(:,nev: dim_xp-1)
    void **CG_p = workspace->CG_p, **evec = workspace->evec, **V_tmp = workspace->V_tmp;
    GCGE_INT    i = 0, j = 0;
    //end[0]: 表示需要正交化向量的最终位置
    //先计算 evec = B*V（：，1：nev），可以在正交化的过程中一直存着
    for(j=0;j<nev;j++)
    {
        ops->GetVecFromMultiVec(V, j, &vec_current);
        ops->GetVecFromMultiVec(evec,j, &vec_tmp);
        if(B == NULL)
        {
            ops->VecAxpby(1.0, vec_current, 0.0, vec_tmp);          
        }
        else
        {
            ops->MatDotVec(B, vec_current, vec_tmp);          
        }
        ops->RestoreVecForMultiVec(evec, j, &vec_tmp);
        ops->RestoreVecForMultiVec(V, j, &vec_current);
    }//end for(j=0;j<nev;j++)
    //用V_tmp存储s剩下来的 B*V(:,nev+1:dim_xp)
    for(j=nev;j<dim_xp;j++)
    {
        ops->GetVecFromMultiVec(V, j, &vec_current);
        ops->GetVecFromMultiVec(V_tmp,j-nev, &vec_tmp);
        if(B == NULL)
        {
            ops->VecAxpby(1.0, vec_current, 0.0, vec_tmp);          
        }
        else
        {
            ops->MatDotVec(B, vec_current, vec_tmp);	  	
        }
        ops->RestoreVecForMultiVec(V_tmp, j-nev, &vec_tmp);
        ops->RestoreVecForMultiVec(V, j, &vec_current);
    }//end for(j=nev;j<dim_xp;j++)

    //开始进行正交化
    //默认做两次正交化，为了增加数值稳定性
    for(i=0; i < 2; ++i)
    {   
        w_length = *end - start;
        for(current = 0; current < nev; ++current)
        {
            // 计算 d_tmp = V^T * vec_tmp
            // 即为 d_tmp[j] = (v_j, vec_current)_B
            mv_s[0] = start;
            mv_e[0] = *(end);
            mv_s[1] = current;
            mv_e[1] = current + 1;
            //d_tmp = vec_tmp*V(:,start:end)
            //统一做内积
            //lda = w_length, 表示矩阵的行数，第一个向量组的列数
            ops->MultiVecInnerProd(V,evec, d_tmp, "ns", mv_s, mv_e, w_length, ops);      
            //减去相应的分量
            for(j=start;j<*(end);j++)
            {
                ops->GetVecFromMultiVec(V, current, &vec_current);
                ops->GetVecFromMultiVec(V, j, &vec_tmp);
                //vec_tmp = vec_tmp - d_tmp[j-start]*vec_current
                //printf("j=%d, xty = %e\n",j, d_tmp[j-start]);
                ops->VecAxpby(-d_tmp[j-start], vec_current, 1.0, vec_tmp);
                ops->RestoreVecForMultiVec(V, j, &vec_tmp);
                ops->RestoreVecForMultiVec(V, current, &vec_current);
            }
        }//end for(current = 0; current < nev; ++current)

        for(current = nev; current < dim_xp; ++current)
        {
            // 计算 d_tmp = V^T * vec_tmp
            // 即为 d_tmp[j] = (v_j, vec_current)_B
            mv_s[0] = start;
            mv_e[0] = *(end);
            mv_s[1] = current - nev;
            mv_e[1] = current + 1 - nev;
            //d_tmp = vec_tmp*V(:,start:end)
            //w_length = mv_e[1] - mv_s[1];
            //统一做内积
            ops->MultiVecInnerProd(V, V_tmp, d_tmp, "ns", mv_s, mv_e, w_length, ops);      
            //减去相应的分量
            for(j=start;j<*(end);j++)
            {
                ops->GetVecFromMultiVec(V, current, &vec_current);
                ops->GetVecFromMultiVec(V, j, &vec_tmp);
                //vec_tmp = vec_tmp - d_tmp[j-start] *vec_current 
                ops->VecAxpby(-d_tmp[j-start], vec_current, 1.0, vec_tmp);
                ops->RestoreVecForMultiVec(V, j, &vec_tmp);
                ops->RestoreVecForMultiVec(V, current, &vec_current);
            }
        }//end for for(current = nev; current < dim_xp; ++current)

        //现在对W本身进行正交化 
        for(current=start;current<*(end);current++)
        {
            //vec_current = W(:,current)
            ops->GetVecFromMultiVec(V, current, &vec_current);
            ops->GetVecFromMultiVec(CG_p, 0, &vec_tmp);
            if(B == NULL)
            {
                ops->VecAxpby(1.0, vec_current, 0.0, vec_tmp);	    
            }
            else
            {
                ops->MatDotVec(B, vec_current, vec_tmp);	  		    
            }
            ops->RestoreVecForMultiVec(CG_p, 0, &vec_tmp);
            ops->RestoreVecForMultiVec(V, current, &vec_current);
            // 计算 d_tmp = V^T * vec_tmp
            // 即为 d_tmp[j] = (v_j, vec_current)_B
            mv_s[0] = current;
            mv_e[0] = *(end);
            mv_s[1] = 0;
            mv_e[1] = 1;
            //d_tmp = vec_tmp*V(:,start:end)
            w_length = mv_e[0] - mv_s[0];
            //统一做内积
            ops->MultiVecInnerProd(V, CG_p, d_tmp, "ns", mv_s, mv_e, w_length, ops);    
            //V[current]的向量范数
            norm_value = sqrt(d_tmp[0]);	  
            //下面看看目前处理的向量是否是一个零向量    	  
            if(norm_value > (orth_para->orth_zero_tol))
            {
                //如果vec_current不为0，就做归一化
                ops->GetVecFromMultiVec(V, current, &vec_current);
                //vec_current = 1/norm_value * vec_current
                ops->VecAxpby(0.0, vec_current, 1.0/norm_value, vec_current);	
                ops->RestoreVecForMultiVec(V, current, &vec_current);
                for(j=current+1;j<*(end);j++)
                {
                    ops->GetVecFromMultiVec(V, current, &vec_current);
                    ops->GetVecFromMultiVec(V, j, &vec_tmp);
                    //vec_tmp = vec_tmp - d_tmp[j-current]/norm_value *vec_current
                    ops->VecAxpby(-d_tmp[j-current]/norm_value, vec_current, 1.0, vec_tmp);
                    ops->RestoreVecForMultiVec(V, j, &vec_tmp);
                    ops->RestoreVecForMultiVec(V, current, &vec_current);
                }//end for(j=current+1;j<*(end);j++)
            }
            else
            {
                //如果vec_current为0，就把此向量与最后一个向量进行指针交换
                //同时(*end)--,把总向量的个数减1
                //同时current--,把当前向量的编号减1,循环到下次仍计算当前编号的向量(即原end)
                mv_s[0] = current; //表示swap中第一个向量组的起始位置
                mv_s[1] = (*end)-1; //表示swap中第二个向量组的起始位置
                mv_e[0] = current + 1; //表示swap中第一个向量组的终止位置
                mv_e[1] = *end; //表示swap中第二个向量组的终止位置
                //向量移动（由用户提供），目的: V(:,current:end-1) = V(:, )
                ops->MultiVecSwap(V, V, mv_s, mv_e, ops);
                (*end)--;
                current--;
                if(orth_para->print_orth_zero == 1)
                {
                    GCGE_Printf("In Orthonormalization, there is a zero vector!, "
                            "current = %d, start = %d, end = %d\n", current, start, *end);		      
                }//end if(orth_para->print_orth_zero == 1)	    
            }//end if(vout > orth_para->orth_zero_tol)	  
        }//end for(current=dim_xp;current<*(end);current++)      
    }//end for(i=0; i <2; ++i)
}//end for the block orthogonalization

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
       }while(ratio < orth_para->reorth_tol && reorth_count < orth_para->max_reorth_time 
 *     && vout > orth_para->orth_zero_tol);
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
void GCGE_OrthonormalizationInSubspace(double *V, GCGE_INT ldV, GCGE_INT nrows, GCGE_INT start, GCGE_INT *end, 
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
        vout = GCGE_ArrayNormInSubspace(V+current*ldV, nrows);
        if(current != 0)
        {
            reorth_count = 0;
            do{
                vin = vout;
                for(idx = 0; idx < current; idx++)
                {
                    //修正的Gram-Schmidt正交化方法
                    //计算ip = (v_idx, v_current)
                    ip = GCGE_ArrayDotArrayInSubspace(V+idx*ldV, V+current*ldV, nrows);
                    //计算v_current = v_current - ip * v_idx
                    GCGE_ArrayAXPBYInSubspace(-ip, V+idx*ldV, 1.0, V+current*ldV, nrows);
                }
                // 计算当前向量的范数
                vout = GCGE_ArrayNormInSubspace(V+current*ldV, nrows);

                ratio = vout/vin;

                reorth_count++;

            }while((ratio < orth_para->reorth_tol) && (reorth_count < orth_para->max_reorth_time)
	            && (vout > orth_para->orth_zero_tol));
        }

        if(vout > orth_para->orth_zero_tol)
        {
            //如果v_current不为0，就做归一化
            GCGE_ArrayScaleInSubspace(1.0/vout, V+current*ldV, nrows);
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
                GCGE_ArrayCopyInSubspace(V+(*end -1)*ldV, V+current*ldV, nrows);
            }//end if (current != *end-1)
            (*end)--;
            current--;
            //如果这个时候需要输出提醒信息，就输出
            if(orth_para->print_orth_zero == 1)
            {
                GCGE_Printf("In OrthonormalizationInSubspace, there is a zero vector!, "
                        "current = %d, start = %d, end = %d\n", current, start, *end);
            }//end if (orth_para->print_orth_zero == 1)
        }//end if (vout > orth_para->orth_zero_tol)
    }//end for current to the end
}//end for this subprogram
#if 0
void GCGE_OrthonormalizationInSubspace(double *V, GCGE_INT ldV, GCGE_INT start, GCGE_INT *end, 
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
        vout = GCGE_ArrayNormInSubspace(V+current*ldV, ldV);
        if(current != 0)
        {
            reorth_count = 0;
            do{
                vin = vout;
                for(idx = 0; idx < current; idx++)
                {
                    //修正的Gram-Schmidt正交化方法
                    //计算ip = (v_idx, v_current)
                    ip = GCGE_ArrayDotArrayInSubspace(V+idx*ldV, V+current*ldV, ldV);
                    //计算v_current = v_current - ip * v_idx
                    GCGE_ArrayAXPBYInSubspace(-ip, V+idx*ldV, 1.0, V+current*ldV, ldV);
                }
                // 计算当前向量的范数
                vout = GCGE_ArrayNormInSubspace(V+current*ldV, ldV);

                ratio = vout/vin;

                reorth_count++;

            }while((ratio < orth_para->reorth_tol) && (reorth_count < orth_para->max_reorth_time)
	            && (vout > orth_para->orth_zero_tol));
        }

        if(vout > orth_para->orth_zero_tol)
        {
            //如果v_current不为0，就做归一化
            GCGE_ArrayScaleInSubspace(1.0/vout, V+current*ldV, ldV);
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
                GCGE_ArrayCopyInSubspace(V+(*end -1)*ldV, V+current*ldV, ldV);
            }//end if (current != *end-1)
            (*end)--;
            current--;
            //如果这个时候需要输出提醒信息，就输出
            if(orth_para->print_orth_zero == 1)
            {
                printf("In OrthonormalizationInSubspace, there is a zero vector!, "
                        "current = %d, start = %d, end = %d\n", current, start, *end);
            }//end if (orth_para->print_orth_zero == 1)
        }//end if (vout > orth_para->orth_zero_tol)
    }//end for current to the end
}//end for this subprogram
#endif

//对V的某些列进行正交化，start: 表示开始的位置，end[0]表示需要正交化向量的最终位置(也返回这个最终位置)，
//矩阵B表示进行内积计算的矩阵，orth_para: 定义正交化方式的参数
//evec, V_tmp: 用来存储正交化过程中 [B*X, B*P] 的临时向量, 
//d_tmp: 用来存储待正交化的向量与之前已经正交化的一个向量的内积 
/** 
 *  V:         input,        要正交化的向量组
 *  start:     input,        要正交化的向量组在V中的起始位置
 *  end:       input|output, 要正交化的向量组在V中的终止位置
 *  B:         input,        正交化矩阵, 可以是NULL
 *  ops:       input,        操作
 *  para:      input,        参数
 *  workspace: input,        工作空间
 *
 *  块正交化过程：
 *  
 *  V = [V1, V2, V3], V中有三个部分, 其中V1,V2部分已经正交, 
 *  V1中有 size_V1(这里是nev) 个向量, V2中有 size_V2(这里是dim_xp-nev) 个向量
 *  V3中有 size_V3(这里是w_length=end-start) 个向量
 *
 *  临时存储空间使用workspace中的 evec(nev个), V_tmp(dim_xp-nev个), CG_p(1个)
 *
 *  下面的过程进行两次(相当于做一次重正交化)
 *  1. 去掉 V3 中 V1,V2 方向的分量
 *     V3 = V3 - [V1,V2] * ([V1,V2]^T * B * V3)
 *   > 计算 V_tmp = B * V3 
 *          V_tmp: w_length 列
 *   > 计算 subspace_dtmp = [V1,V2]^T * V_tmp
 *          subspace_dtmp: GCGE_DOUBLE *数组
 *            size_V1+size_V2 行, 即dim_xp行,w_length 列
 *   > 计算 V_tmp = [V1,V2] * subspace_dtmp
 *          V_tmp: w_length 列
 *   > 计算 V3 = V3 - V_tmp
 *
 *  2. V3 自身做正交化 
 *   > 计算 M = V3^T * B *V3
 *   > 计算 M C = C Theta
 *   > 计算 y_i = 1/sqrt(Theta_i) * y_i
 *   > 计算 V3 = V3 * Y
 *
 */
//Classical Block Orthonormalizationization
void GCGE_CBOrthonormalization(void **V, GCGE_INT start, GCGE_INT *end, 
                      void *B, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
    GCGE_ORTH_PARA *orth_para = para->orth_para; 
    GCGE_INT        current = 0; //表示for循环中当前正在进行正交化的向量编号
    GCGE_INT        w_length = (*end) - start; //表示W向量的个数
    GCGE_INT        reorth_count = 0; //统计(重)正交化次数
    GCGE_INT        mv_s[2] = { 0, 0 }; //使用向量组的初始位置
    GCGE_INT        mv_e[2] = { 0, 0 }; //使用向量组的终止位置

    GCGE_INT        dim_x = workspace->dim_x;
    GCGE_INT        dim_xp = workspace->dim_xp;
    GCGE_INT        nev = para->nev;

    GCGE_DOUBLE     norm_value = 0.0;  //存储当前向量的范数
    GCGE_DOUBLE    *d_tmp = workspace->subspace_dtmp;
    void           *vec_current; //当前正在进行正交化的向量
    void           *vec_tmp;
    void          **CG_p = workspace->CG_p;
    void          **evec = workspace->evec;
    void          **V_tmp = workspace->V_tmp;
    GCGE_INT        i = 0, j = 0;

    //求解特征值问题所用参数
    GCGE_DOUBLE     vl = 0.0;
    GCGE_DOUBLE     vu = 0.0;
    GCGE_DOUBLE     abstol = 0.0;
    GCGE_DOUBLE    *tmp_evec = workspace->subspace_evec;
    GCGE_DOUBLE    *tmp_eval = tmp_evec + w_length*w_length;
    GCGE_DOUBLE    *dwork_space = d_tmp + w_length*w_length;
    GCGE_INT        il = 1;
    GCGE_INT        iu = w_length;
    GCGE_INT        m = iu-il+1;
    GCGE_INT        lwork = 8*w_length;
    GCGE_INT        liwork = 5*w_length;
    GCGE_INT       *subspace_itmp = workspace->subspace_itmp;
    GCGE_INT       *isuppz = subspace_itmp + 10*w_length; 
    GCGE_INT       *ifail = subspace_itmp + 5*w_length; 
    GCGE_INT        info = 0;

    GCGE_INT        tmp_eval_nonzero_start = 0; //用于做线性组合的小规模特征向量的起始位置
    GCGE_INT        nonzero_flag = 0;//用于确定内积是否都接近0

    //开始进行正交化
    //默认做两次正交化，为了增加数值稳定性
    for(i=0; i < 2; ++i)
    {   
        w_length = *end - start;

        //计算 V_tmp = B * V3 
        if(B == NULL)
        {
            mv_s[0] = start;
            mv_e[0] = *end;
            mv_s[1] = 0;
            mv_e[1] = w_length;
            ops->MultiVecAxpby(1.0, V, 0.0, V_tmp, mv_s, mv_e, ops);	  	
        }
        else
        {
            mv_s[0] = start;
            mv_e[0] = *end;
            mv_s[1] = 0;
            mv_e[1] = w_length;
            ops->MatDotMultiVec(B, V, V_tmp, mv_s, mv_e, ops);
        }

        //计算 subspace_dtmp = [V1,V2]^T * V_tmp
        mv_s[0] = 0;
        mv_e[0] = start;
        mv_s[1] = 0;
        mv_e[1] = w_length;
        ops->MultiVecInnerProd(V, V_tmp, d_tmp, "ns", mv_s, mv_e, start, ops);    

        //计算 V_tmp = [V1,V2] * subspace_dtmp
        mv_s[0] = 0;
        mv_e[0] = start;
        mv_s[1] = 0;
        mv_e[1] = w_length;
        ops->MultiVecLinearComb(V, V_tmp, mv_s, mv_e, d_tmp, start, NULL, -1, ops);

        //计算 V3 = V3 - V_tmp
        mv_s[0] = 0;
        mv_e[0] = w_length;
        mv_s[1] = start;
        mv_e[1] = *end;
        ops->MultiVecAxpby(-1.0, V_tmp, 1.0, V, mv_s, mv_e, ops);
    }

    //现在对W本身进行正交化,因为这种正交化方式只适用于精度要求低的情况
    //所以这里暂时不做重正交化
    GCGE_INT rank = 0;
    for(i=0; i < 1; ++i)
    {   
        //计算 M = V3^T * B *V3
        //计算 V_tmp = B * V3 
        for(j=0;j<w_length;j++)
        {
            ops->GetVecFromMultiVec(V, j+start, &vec_current);
            ops->GetVecFromMultiVec(V_tmp, j, &vec_tmp);
            if(B == NULL)
            {
                ops->VecAxpby(1.0, vec_current, 0.0, vec_tmp);	  	
            }
            else
            {
                ops->MatDotVec(B, vec_current, vec_tmp);	  	
            }
            ops->RestoreVecForMultiVec(V_tmp, j, &vec_tmp);
            ops->RestoreVecForMultiVec(V, j+start, &vec_current);
        }//end for(j=0;j<w_length;j++)

        //计算 subspace_dtmp = V3^T * V_tmp
        mv_s[0] = start;
        mv_e[0] = *end;
        mv_s[1] = 0;
        mv_e[1] = w_length;
        ops->MultiVecInnerProd(V, V_tmp, d_tmp, "ns", mv_s, mv_e, w_length, ops);    

        nonzero_flag = 0;
        for(j=0; j<w_length*w_length; j++)
        {
            if(d_tmp[j] > orth_para->orth_zero_tol)
            {
                nonzero_flag = 1;
                break;
            }
        }
        if(nonzero_flag)
        {
            //计算 M C = C Theta
            iu = w_length;
            m = iu - il + 1;
            if(para->opt_bcast)
            {
                if(para->use_mpi_bcast)
                {
#if GCGE_USE_MPI
                    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
                }
            }
            //此时, 如果使用 MPI 且要用 MPI_Bcast, 那么 非0号 进程 rank != 0, 
            //此时 0号 进程 rank==0, 如果 不用MPI_Bcast 也是 rank != 0
            //那么, 如果 rank==0, 就计算特征值问题, 否则不用计算, 等待广播
            if(rank == 0)
            {
                ops->DenseMatEigenSolver("V", "I", "U", &w_length, d_tmp, 
                        &w_length, &vl, &vu, &il, &iu, &abstol, 
                        &m, tmp_eval, tmp_evec, &w_length, 
                        isuppz, dwork_space, &lwork,
                        subspace_itmp, &liwork, ifail, &info);
            }
            if(para->use_mpi_bcast)
            {
#if GCGE_USE_MPI
                memcpy(tmp_evec+iu*w_length, tmp_eval, m*sizeof(GCGE_DOUBLE));
                MPI_Bcast(tmp_evec, m*w_length+m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                memcpy(tmp_eval, tmp_evec+iu*w_length, m*sizeof(GCGE_DOUBLE));
#endif
            }

            //计算 y_i = 1/sqrt(Theta_i) * y_i
            tmp_eval_nonzero_start = 0;
            for(j=0; j<w_length; j++)
            {
                if(tmp_eval[j] < orth_para->orth_zero_tol)
                {
                    //表明这个位置是0，那非0的要从下一个开始
                    tmp_eval_nonzero_start = j+1;
                }
                else
                {
                    GCGE_ArrayScaleInSubspace(1.0/sqrt(tmp_eval[j]), tmp_evec+j*w_length, w_length);
                }
            }
            if(tmp_eval_nonzero_start)
            {
                GCGE_Printf("tmp_eval_nonzero_start: %d\n", tmp_eval_nonzero_start);
                for(j=0; j<w_length; j++)
                    GCGE_Printf("tmp_eval[%d]: %e\n", j, tmp_eval[j]);
            }

            //w_length为新的非零向量个数
            w_length = w_length - tmp_eval_nonzero_start;

            //计算 V3 = V3 * Y
            //计算 V_tmp = [V1,V2] * subspace_dtmp
            mv_s[0] = start;
            mv_e[0] = *end;
            mv_s[1] = 0;
            mv_e[1] = w_length;
            ops->MultiVecLinearComb(V, V_tmp, mv_s, mv_e, 
                    tmp_evec+tmp_eval_nonzero_start*w_length, w_length, NULL, -1, ops);
            *end = start + w_length;
            //计算 V3 = V_tmp
            mv_s[0] = start;
            mv_e[0] = *end;
            mv_s[1] = 0;
            mv_e[1] = w_length;
            ops->MultiVecSwap(V, V_tmp, mv_s, mv_e, ops);
        }
        else
        {
            GCGE_Printf("W all zero!\n");
            *end = start;
            i = 10;
        }
    }//end for(i=0; i <2; ++i)
}//end for the block orthogonalization

//Classical Block Orthonormalizationization
void GCGE_SCBOrthonormalization(void **V, GCGE_INT start, GCGE_INT *end, 
                      void *B, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
    GCGE_ORTH_PARA  *orth_para = para->orth_para; 
    GCGE_INT        current = 0;               //表示for循环中当前正在进行正交化的向量编号
    GCGE_INT        w_length = (*end) - start; //表示W向量的个数
    GCGE_INT        reorth_count = 0; //统计(重)正交化次数
    GCGE_INT        mv_s[2] = { 0, 0 }; //使用向量组的初始位置
    GCGE_INT        mv_e[2] = { 0, 0 }; //使用向量组的终止位置

    GCGE_INT        dim_x  = workspace->dim_x;
    GCGE_INT        dim_xp = workspace->dim_xp;
    GCGE_INT        nev = para->nev;

    GCGE_DOUBLE     norm_value = 0.0;  //存储当前向量的范数
    GCGE_DOUBLE    *d_tmp = workspace->subspace_dtmp;
    void           *vec_current; //当前正在进行正交化的向量
    void           *vec_tmp;
    void          **CG_p = workspace->CG_p;
    //void          **evec = workspace->evec;
    void           **V_tmp = workspace->V_tmp;
    GCGE_INT        i = 0, j = 0;

    //求解特征值问题所用参数
    GCGE_DOUBLE     vl = 0.0;
    GCGE_DOUBLE     vu = 0.0;
    GCGE_DOUBLE     abstol = 0.0;
    GCGE_DOUBLE    *tmp_evec = workspace->subspace_evec;
    //tmp_evec前 dim_W*dim_W的位置用来存储特真向量
    GCGE_DOUBLE    *tmp_eval = tmp_evec + w_length*w_length;
    GCGE_DOUBLE    *dwork_space = d_tmp + w_length*w_length;
    GCGE_INT        il = 1;
    GCGE_INT        iu = w_length;
    GCGE_INT        m = iu-il+1;
    GCGE_INT        lwork = 8*w_length;
    GCGE_INT        liwork = 5*w_length;
    GCGE_INT       *subspace_itmp = workspace->subspace_itmp;
    GCGE_INT       *isuppz = subspace_itmp + 10*w_length; 
    GCGE_INT       *ifail = subspace_itmp + 5*w_length; 
    GCGE_INT        info = 0;
    GCGE_DOUBLE     Orth_Tol = orth_para->scbgs_reorth_tol;
    GCGE_DOUBLE     value_inner = 1.0, tmp;
    GCGE_INT        length = start *w_length, iter = 0;

    GCGE_INT        tmp_eval_nonzero_start = 0; //用于做线性组合的小规模特征向量的起始位置
    GCGE_INT        nonzero_flag = 0;           //用于确定内积是否都接近0

    //开始进行正交化, 首先把W中含有的X和P的分量都减去，这一部分需要进行快操作，这样才有效率
    //默认做两次正交化，为了增加数值稳定性,
    while(value_inner > Orth_Tol)
    {   
        //w_length = *end - start;
        iter ++;      
        
        //计算 V_tmp = B * V3 
        if(B == NULL)
        {
            mv_s[0] = start;
            mv_e[0] = *end;
            mv_s[1] = 0;
            mv_e[1] = w_length;
            ops->MultiVecAxpby(1.0, V, 0.0, V_tmp, mv_s, mv_e, ops);	  	
        }
        else
        {
            mv_s[0] = start;
            mv_e[0] = *end;
            mv_s[1] = 0;
            mv_e[1] = w_length;
            ops->MatDotMultiVec(B, V, V_tmp, mv_s, mv_e, ops);
        }
        //计算 subspace_dtmp = [V1,V2]^T * V_tmp
        mv_s[0] = 0;
        mv_e[0] = start;
        mv_s[1] = 0;
        mv_e[1] = w_length;
        ops->MultiVecInnerProd(V, V_tmp, d_tmp, "ns", mv_s, mv_e, start, ops);  
        value_inner = 0.0;
        for(i=0;i<length;i++)
        { 
            tmp = d_tmp[i];
            value_inner += tmp*tmp;
         }
         value_inner = sqrt(value_inner);
         GCGE_Printf("the %d-th orth, inner value=%e\n",iter, value_inner);    
         

        //计算 V_tmp = [V1,V2] * subspace_dtmp
        mv_s[0] = 0;
        mv_e[0] = start;
        mv_s[1] = 0;
        mv_e[1] = w_length;
        ops->MultiVecLinearComb(V, V_tmp, mv_s, mv_e, d_tmp, start, NULL, -1, ops);

        //计算 V3 = V3 - V_tmp
        mv_s[0] = 0;
        mv_e[0] = w_length;
        mv_s[1] = start;
        mv_e[1] = *end;
        ops->MultiVecAxpby(-1.0, V_tmp, 1.0, V, mv_s, mv_e, ops);

        if(iter >= orth_para->max_reorth_time)
            break;
    }//end while for value_inner
    
    
    //对W自身进行正交化
    GCGE_INT step = ((*end - start) < orth_para->w_orth_block_size)?(*end - start) : orth_para->w_orth_block_size;
    // 对V中[start：end]向量组本身进行正交化   
    GCGE_INT w_start, w_end, old_w_end;    
    
    //第一步对开始部分进行正交化
    w_start = start;        
    w_end = w_start + step;
    if(w_end> *end)
        w_end = *end;
    //接下来就是对 W[w_start:w_end]进行正交化
    GCGE_SCBOrth_Self(V, w_start, &w_end, B, ops, para, workspace);
    
    //然后进行迭代，在每一步首先减去前面的分量，然后再对自己的部分进行正交化
    //current = w_end;
    while(w_end < *end)
    {
        //第一步：把之前的分量减去
        w_start = w_end;
        w_end = w_start + step;
        if(w_end> *end)
             w_end = *end;
        old_w_end = w_end;
        GCGE_SCBOrth_Minus(V, w_start, &w_end, B, ops, para, workspace);
        //第二步：对自身进行正交化: 如果中间出现零向量，会返回新的w_end和workspace->dim_xpw
        GCGE_SCBOrth_Self(V, w_start, &w_end, B, ops, para, workspace);
        //如果出现零向量，把V的最后几列与现在w_end之后的几列相交换，并且同时end下降相应的
        if(w_end<old_w_end)
        {
           *end -= old_w_end - w_end;
           mv_s[0] = w_end;      //表示swap中第一个向量组的起始位置
           mv_s[1] = *end;     //表示swap中第二个向量组的起始位置
           mv_e[0] = old_w_end;  //表示swap中第一个向量组的终止位置
           mv_e[1] = *end + old_w_end - w_end;         //表示swap中第二个向量组的终止位置
           ops->MultiVecSwap(V, V, mv_s, mv_e, ops);
        }//end for if w_end<old_w_end              
    }//end for(current = dim_xpw; current<*(end); current++)      
}//end for the block orthogonalization


//subfunction: 把当前的向量组见去之前已经正交化好的向量组
//注意，这里是一个定制的正交化过程，现在是对V的W部分进行正交化，之前的X和P部分已经处理完备，现在
//这里减去的分量都是W的部分
void GCGE_SCBOrth_Minus(void **V, GCGE_INT start, GCGE_INT *end, 
        void *B, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
    GCGE_ORTH_PARA  *orth_para = para->orth_para; 
    GCGE_INT        current = 0; //表示for循环中当前正在进行正交化的向量编号
    GCGE_INT        w_length = (*end) - start; //表示目前处理的向量个数
    GCGE_INT        reorth_count = 0; //统计(重)正交化次数
    GCGE_INT        mv_s[2] = { 0, 0 }; //使用向量组的初始位置
    GCGE_INT        mv_e[2] = { 0, 0 }; //使用向量组的终止位置
    
    GCGE_INT        dim_xp = workspace->dim_xp;

    //GCGE_DOUBLE     norm_value = 0.0;  //存储当前向量的范数
    GCGE_DOUBLE     *d_tmp = workspace->subspace_dtmp;
    void            **V_tmp = workspace->V_tmp;
    GCGE_INT        i = 0;
    GCGE_DOUBLE     Orth_Tol = orth_para->scbgs_reorth_tol;
    GCGE_DOUBLE     value_inner = 1.0, tmp;
    GCGE_INT        length = (start - dim_xp) *w_length, iter = 0;
    //开始进行正交化
    //默认做两次正交化，为了增加数值稳定性
    while(value_inner > Orth_Tol)
    {   
        //w_length = *end - start;
        iter ++;        
        //计算 V_tmp = B * V3 
        if(B == NULL)
        {
            mv_s[0] = start;
            mv_e[0] = *end;
            mv_s[1] = 0;
            mv_e[1] = w_length;
            ops->MultiVecAxpby(1.0, V, 0.0, V_tmp, mv_s, mv_e, ops);	  	
        }
        else
        {
            mv_s[0] = start;
            mv_e[0] = *end;
            mv_s[1] = 0;
            mv_e[1] = w_length;
            ops->MatDotMultiVec(B, V, V_tmp, mv_s, mv_e, ops);
        }
        //计算 subspace_dtmp = W1^T * V_tmp
        mv_s[0] = dim_xp;
        mv_e[0] = start;
        mv_s[1] = 0;
        mv_e[1] = w_length;
        ops->MultiVecInnerProd(V, V_tmp, d_tmp, "ns", mv_s, mv_e, start, ops);  
        value_inner = 0.0;
        for(i=0;i<length;i++)
        { 
            tmp = d_tmp[i];
            value_inner += tmp*tmp;
         }
         value_inner = sqrt(value_inner);
         GCGE_Printf("the %d-th orth, inner value=%e\n",iter, value_inner);    
         

        //计算 V_tmp = W1 * subspace_dtmp
        mv_s[0] = dim_xp;
        mv_e[0] = start;
        mv_s[1] = 0;
        mv_e[1] = w_length;
        ops->MultiVecLinearComb(V, V_tmp, mv_s, mv_e, d_tmp, start, NULL, -1, ops);

        //计算 W2 = W2 - V_tmp
        mv_s[0] = 0;
        mv_e[0] = w_length;
        mv_s[1] = start;
        mv_e[1] = *end;
        ops->MultiVecAxpby(-1.0, V_tmp, 1.0, V, mv_s, mv_e, ops);

        if(iter >= orth_para->max_reorth_time)
            break;
    }//end while for value_inner 
}//end for orth_minus  
        

//subfunction for orthogonalization process: Orthonormalizationization for the vectors self
//对V中[start：end]向量组本身进行正交化    
void GCGE_SCBOrth_Self(void **V, GCGE_INT start, GCGE_INT *end, 
              void *B, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
   GCGE_ORTH_PARA  *orth_para = para->orth_para; 
   GCGE_INT        iter = 0, current, w_lengths;
   GCGE_INT        w_length;
   GCGE_INT        j = 0;
   GCGE_DOUBLE     norm_value = 0.0;  //存储当前向量的范数
   GCGE_DOUBLE     *d_tmp = workspace->subspace_dtmp;
   void            *vec_current; //当前正在进行正交化的向量
   void            *vec_tmp;
   void            **CG_p = workspace->CG_p;
   GCGE_INT        mv_s[2] = { 0, 0 }; //使用向量组的初始位置
   GCGE_INT        mv_e[2] = { 0, 0 }; //使用向量组的终止位置
   //现在对W本身进行正交化 
   for(iter=0;iter<(orth_para->scbgs_wself_max_reorth_time);iter++)
   {
      for(current=start;current<*(end);current++)
      {
           //vec_current = W(:,current)
           ops->GetVecFromMultiVec(V, current, &vec_current);
           ops->GetVecFromMultiVec(CG_p, 0, &vec_tmp);
           if(B == NULL)
           {
               ops->VecAxpby(1.0, vec_current, 0.0, vec_tmp);	    
           }
           else
           {
               ops->MatDotVec(B, vec_current, vec_tmp);	  		    
           }
           ops->RestoreVecForMultiVec(CG_p, 0, &vec_tmp);
           ops->RestoreVecForMultiVec(V, current, &vec_current);
           // 计算 d_tmp = V^T * vec_tmp
           // 即为 d_tmp[j] = (v_j, vec_current)_B
           mv_s[0] = current;
           mv_e[0] = *(end);
           mv_s[1] = 0;
           mv_e[1] = 1;
           //d_tmp = vec_tmp*V(:,start:end)
           w_length = mv_e[0] - mv_s[0];
           //统一做内积
           ops->MultiVecInnerProd(V, CG_p, d_tmp, "ns", mv_s, mv_e, w_length, ops);    
           //V[current]的向量范数
           norm_value = sqrt(d_tmp[0]);	  
           //下面看看目前处理的向量是否是一个零向量    	  
           if(norm_value > (orth_para->orth_zero_tol))
           {
               //如果vec_current不为0，就做归一化
               ops->GetVecFromMultiVec(V, current, &vec_current);
               //vec_current = 1/norm_value * vec_current
               ops->VecAxpby(0.0, vec_current, 1.0/norm_value, vec_current);	
               ops->RestoreVecForMultiVec(V, current, &vec_current);
               for(j=current+1;j<*(end);j++)
               {
                    ops->GetVecFromMultiVec(V, current, &vec_current);
                    ops->GetVecFromMultiVec(V, j, &vec_tmp);
                    //vec_tmp = vec_tmp - d_tmp[j-current]/norm_value *vec_current
                    ops->VecAxpby(-d_tmp[j-current]/norm_value, vec_current, 1.0, vec_tmp);
                    ops->RestoreVecForMultiVec(V, j, &vec_tmp);
                    ops->RestoreVecForMultiVec(V, current, &vec_current);
                }//end for(j=current+1;j<*(end);j++)
            }
            else
            {
                //如果vec_current为0，就把此向量与最后一个向量进行指针交换
                //同时(*end)--,把总向量的个数减1
                //同时current--,把当前向量的编号减1,循环到下次仍计算当前编号的向量(即原end)
                mv_s[0] = current;      //表示swap中第一个向量组的起始位置
                mv_s[1] = (*end)-1;     //表示swap中第二个向量组的起始位置
                mv_e[0] = current + 1;  //表示swap中第一个向量组的终止位置
                mv_e[1] = *end;         //表示swap中第二个向量组的终止位置
                //向量移动(由用户提供), 目的: V(:,current:end-1) = V(:, )
                ops->MultiVecSwap(V, V, mv_s, mv_e, ops);
                (*end)--;
                current--;                          
                if(orth_para->print_orth_zero == 1)
                {
                    GCGE_Printf("In Orthonormalization, there is a zero vector!, "
                            "current = %d, start = %d, end = %d\n", current, start, *end);		      
                }//end if(orth_para->print_orth_zero == 1)	    
            }//end if(vout > orth_para->orth_zero_tol)	  
        }//end for(current=dim_xp;current<*(end);current++)      
    }//end for(i=0; i <2; ++i)
}//GCGE_SCBOrth_Self


#if 0
//Classical Block Orthonormalizationization
void GCGE_SCBOrthonormalization(void **V, GCGE_INT start, GCGE_INT *end, 
                      void *B, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
    GCGE_ORTH_PARA  *orth_para = para->orth_para; 
    GCGE_INT        current = 0; //表示for循环中当前正在进行正交化的向量编号
    GCGE_INT        w_length = (*end) - start; //表示W向量的个数
    GCGE_INT        reorth_count = 0; //统计(重)正交化次数
    GCGE_INT        mv_s[2] = { 0, 0 }; //使用向量组的初始位置
    GCGE_INT        mv_e[2] = { 0, 0 }; //使用向量组的终止位置

    GCGE_INT        dim_x  = workspace->dim_x;
    GCGE_INT        dim_xp = workspace->dim_xp;
    GCGE_INT        nev = para->nev;

    GCGE_DOUBLE     norm_value = 0.0;  //存储当前向量的范数
    GCGE_DOUBLE    *d_tmp = workspace->subspace_dtmp;
    void           *vec_current; //当前正在进行正交化的向量
    void           *vec_tmp;
    void          **CG_p = workspace->CG_p;
    //void          **evec = workspace->evec;
    void           **V_tmp = workspace->V_tmp;
    GCGE_INT        i = 0, j = 0;

    //求解特征值问题所用参数
    GCGE_DOUBLE     vl = 0.0;
    GCGE_DOUBLE     vu = 0.0;
    GCGE_DOUBLE     abstol = 0.0;
    GCGE_DOUBLE    *tmp_evec = workspace->subspace_evec;
    //tmp_evec前 dim_W*dim_W的位置用来存储特真向量
    GCGE_DOUBLE    *tmp_eval = tmp_evec + w_length*w_length;
    GCGE_DOUBLE    *dwork_space = d_tmp + w_length*w_length;
    GCGE_INT        il = 1;
    GCGE_INT        iu = w_length;
    GCGE_INT        m = iu-il+1;
    GCGE_INT        lwork = 8*w_length;
    GCGE_INT        liwork = 5*w_length;
    GCGE_INT       *subspace_itmp = workspace->subspace_itmp;
    GCGE_INT       *isuppz = subspace_itmp + 10*w_length; 
    GCGE_INT       *ifail = subspace_itmp + 5*w_length; 
    GCGE_INT        info = 0;
    GCGE_DOUBLE     Orth_Tol = orth_para->scbgs_reorth_tol;
    GCGE_DOUBLE     value_inner = 1.0, tmp;
    GCGE_INT        length = start *w_length, iter = 0;

    GCGE_INT        tmp_eval_nonzero_start = 0; //用于做线性组合的小规模特征向量的起始位置
    GCGE_INT        nonzero_flag = 0;//用于确定内积是否都接近0

    //开始进行正交化
    //默认做两次正交化，为了增加数值稳定性
    while(value_inner > Orth_Tol)
    {   
        //w_length = *end - start;
        iter ++;
        
        
        //计算 V_tmp = B * V3 
        if(B == NULL)
        {
            mv_s[0] = start;
            mv_e[0] = *end;
            mv_s[1] = 0;
            mv_e[1] = w_length;
            ops->MultiVecAxpby(1.0, V, 0.0, V_tmp, mv_s, mv_e, ops);	  	
        }
        else
        {
            mv_s[0] = start;
            mv_e[0] = *end;
            mv_s[1] = 0;
            mv_e[1] = w_length;
            ops->MatDotMultiVec(B, V, V_tmp, mv_s, mv_e, ops);
        }

        //计算 subspace_dtmp = [V1,V2]^T * V_tmp
        mv_s[0] = 0;
        mv_e[0] = start;
        mv_s[1] = 0;
        mv_e[1] = w_length;
        ops->MultiVecInnerProd(V, V_tmp, d_tmp, "ns", mv_s, mv_e, start, ops);  
        value_inner = 0.0;
        for(i=0;i<length;i++)
        { 
            tmp = d_tmp[i];
            value_inner += tmp*tmp;
         }
         value_inner = sqrt(value_inner);
         GCGE_Printf("the %d-th orth, inner value=%e\n",iter, value_inner);    
         

        //计算 V_tmp = [V1,V2] * subspace_dtmp
        mv_s[0] = 0;
        mv_e[0] = start;
        mv_s[1] = 0;
        mv_e[1] = w_length;
        ops->MultiVecLinearComb(V, V_tmp, mv_s, mv_e, d_tmp, start, NULL, -1, ops);

        //计算 V3 = V3 - V_tmp
        mv_s[0] = 0;
        mv_e[0] = w_length;
        mv_s[1] = start;
        mv_e[1] = *end;
        ops->MultiVecAxpby(-1.0, V_tmp, 1.0, V, mv_s, mv_e, ops);

        if(iter >= orth_para->max_reorth_time)
            break;
    }//end while for value_inner

    //现在对W本身进行正交化 
    for(iter=0;iter<(orth_para->scbgs_wself_max_reorth_time);iter++)
    {
        for(current=start;current<*(end);current++)
        {
            //vec_current = W(:,current)
            ops->GetVecFromMultiVec(V, current, &vec_current);
            ops->GetVecFromMultiVec(CG_p, 0, &vec_tmp);
            if(B == NULL)
            {
                ops->VecAxpby(1.0, vec_current, 0.0, vec_tmp);	    
            }
            else
            {
                ops->MatDotVec(B, vec_current, vec_tmp);	  		    
            }
            ops->RestoreVecForMultiVec(CG_p, 0, &vec_tmp);
            ops->RestoreVecForMultiVec(V, current, &vec_current);
            // 计算 d_tmp = V^T * vec_tmp
            // 即为 d_tmp[j] = (v_j, vec_current)_B
            mv_s[0] = current;
            mv_e[0] = *(end);
            mv_s[1] = 0;
            mv_e[1] = 1;
            //d_tmp = vec_tmp*V(:,start:end)
            w_length = mv_e[0] - mv_s[0];
            //统一做内积
            ops->MultiVecInnerProd(V, CG_p, d_tmp, "ns", mv_s, mv_e, w_length, ops);    
            //V[current]的向量范数
            norm_value = sqrt(d_tmp[0]);	  
            //下面看看目前处理的向量是否是一个零向量    	  
            if(norm_value > (orth_para->orth_zero_tol))
            {
                //如果vec_current不为0，就做归一化
                ops->GetVecFromMultiVec(V, current, &vec_current);
                //vec_current = 1/norm_value * vec_current
                ops->VecAxpby(0.0, vec_current, 1.0/norm_value, vec_current);	
                ops->RestoreVecForMultiVec(V, current, &vec_current);
                for(j=current+1;j<*(end);j++)
                {
                    ops->GetVecFromMultiVec(V, current, &vec_current);
                    ops->GetVecFromMultiVec(V, j, &vec_tmp);
                    //vec_tmp = vec_tmp - d_tmp[j-current]/norm_value *vec_current
                    ops->VecAxpby(-d_tmp[j-current]/norm_value, vec_current, 1.0, vec_tmp);
                    ops->RestoreVecForMultiVec(V, j, &vec_tmp);
                    ops->RestoreVecForMultiVec(V, current, &vec_current);
                }//end for(j=current+1;j<*(end);j++)
            }
            else
            {
                //如果vec_current为0，就把此向量与最后一个向量进行指针交换
                //同时(*end)--,把总向量的个数减1
                //同时current--,把当前向量的编号减1,循环到下次仍计算当前编号的向量(即原end)
                mv_s[0] = current; //表示swap中第一个向量组的起始位置
                mv_s[1] = (*end)-1; //表示swap中第二个向量组的起始位置
                mv_e[0] = current + 1; //表示swap中第一个向量组的终止位置
                mv_e[1] = *end; //表示swap中第二个向量组的终止位置
                //向量移动（由用户提供），目的: V(:,current:end-1) = V(:, )
                ops->MultiVecSwap(V, V, mv_s, mv_e, ops);
                (*end)--;
                current--;
                if(orth_para->print_orth_zero == 1)
                {
                    GCGE_Printf("In Orthonormalization, there is a zero vector!, "
                            "current = %d, start = %d, end = %d\n", current, start, *end);		      
                }//end if(orth_para->print_orth_zero == 1)	    
            }//end if(vout > orth_para->orth_zero_tol)	  
        }//end for(current=dim_xp;current<*(end);current++)      
    }//end for(i=0; i <2; ++i)
}//end for the block orthogonalization
#endif

/** 
 *  V:         input,        要正交化的向量组
 *  ldV:       input,        子空间向量组V的长度
 *  start:     input,        要正交化的向量组在V中的起始位置
 *  end:       input|output, 要正交化的向量组在V中的终止位置
 *  B:         input,        正交化矩阵, 可以是NULL,这里只考虑B==NULL的情况
 *  ldB:       input,        子空间矩阵B的规模
 *  orth_para: input,        正交化参数
 *  d_tmp:     input,        临时存储空间
 *
 *  块正交化过程：
 *  
 *  V = [V1, V2], V中有两个部分, 其中V1部分已经正交, 
 *  V1中有 size_V1(这里是start)个向量, V2中有 size_V2(这里是end-start) 个向量
 *
 *  临时存储空间使用workspace中的 evec(nev个), V_tmp(dim_xp-nev个), CG_p(1个)
 *
 *  下面的过程进行两次(相当于做一次重正交化)
 *  1. 去掉 V2 中 V1 方向的分量
 *     for current = start_V1:end_V1
 *         计算 d_tmp = (V2, V[current])
 *         for V2_current = start_V2:end_V2
 *             计算 V[V2_current] = V[V2_current] - d_tmp[V2_current-start_V2] * V[current]
 *         end
 *     end
 *
 *  2. V2 自身做正交化 
 *     for current = start_V2:end_V2
 *         计算 d_tmp = (V[current:end_V2], V[current])
 *         norm = sqrt(d_tmp[0])
 *         如果norm > orth_zero_tol
 *            for V2_current = current:end_V2
 *                计算 V[V2_current] = V[V2_current] - d_tmp[V2_current-i]/norm * V[current]
 *            end
 *         否则，判定向量V[V2_current]为0
 *     end
 *
 */
void GCGE_BOrthonormalizationInSubspace(double *V, GCGE_INT ldV, GCGE_INT nrows, GCGE_INT start, 
        GCGE_INT *end, void *B, GCGE_INT ldB, GCGE_ORTH_PARA *orth_para, 
        GCGE_DOUBLE *d_tmp, GCGE_OPS *ops)
{
    GCGE_Printf("Use GCGE_BOrthonormalizationInSubspace.\n");
    GCGE_INT    current = 0;        //当前进行正交化操作的向量编号
    GCGE_INT    current_V2 = 0;     //当前进行正交化操作的向量编号
    GCGE_INT    reorth_count = 0;
    GCGE_INT    orth_length = *end - start;
    GCGE_INT    orth_ncols = 1;
    GCGE_INT    max_reorth_count = orth_para->max_reorth_time;
    GCGE_DOUBLE orth_zero_tol = orth_para->orth_zero_tol;
    GCGE_DOUBLE alpha = 1.0;
    GCGE_DOUBLE beta = 0.0;
    GCGE_DOUBLE norm = 0.0;
    max_reorth_count = 2;
    for(reorth_count=0; reorth_count < max_reorth_count; reorth_count++)
    {
        orth_length = *end - start;
        //去掉 V2 中 V1 方向的分量
        for(current = 0; current < start; ++current)
        {
            //计算 d_tmp = (V2, V[current])
            ops->DenseMatDotDenseMat("T", "N", &orth_length, &orth_ncols,
                    &nrows, &alpha, V+start*ldV, &ldV, V+current*ldV, &ldV,
                    &beta, d_tmp, &orth_length);
            //ops->DenseMatDotVec("T", &ldV, &orth_length, &alpha, V+start*ldV,
            //        &ldV, V+current*ldV, &orth_ncols, &beta, d_tmp, &orth_ncols);
            for(current_V2 = start; current_V2 < *end; current_V2++)
            {
                if(fabs(d_tmp[current_V2-start]) > 10.0*orth_zero_tol)
                     GCGE_ArrayAXPBYInSubspace(-d_tmp[current_V2-start], V+current*ldV,
                        1.0, V+current_V2*ldV, nrows);
            }
        }
        //V2 自身做正交化 
        for(current = start; current < *end; ++current)
        {
            orth_length = *end - current;
            //计算 d_tmp = (V[current:end_V2], V[current])
            ops->DenseMatDotDenseMat("T", "N", &orth_length, &orth_ncols,
                    &nrows, &alpha, V+current*ldV, &ldV, V+current*ldV, &ldV,
                    &beta, d_tmp, &orth_length);
            //ops->DenseMatDotVec("T", &ldV, &orth_length, &alpha, 
            //        V+current*ldV, &ldV, V+current*ldV, &orth_ncols, 
            //        &beta, d_tmp, &orth_ncols);
            norm = sqrt(d_tmp[0]);
            if(norm > orth_zero_tol)
            {
                GCGE_ArrayScaleInSubspace(1.0/norm, V+current*ldV, nrows);
                for(current_V2=current+1; current_V2 < *end; current_V2++)
                {
                    if(fabs(d_tmp[current_V2 - current]/norm) > 10.0*orth_zero_tol)
                        GCGE_ArrayAXPBYInSubspace(-d_tmp[current_V2 - current]/norm, 
                            V+current*ldV, 1.0, V+current_V2*ldV, nrows);
                }
            }
            else
            {
                //如果v_current为0，就把最后一个向量拷贝到当前向量的位置
                //如果本向量就是最后一个向量，那就不拷贝(考虑如果是最后一个直接跳出)
                //同时(*end)--,把总向量的个数减1
                //同时current--,把当前向量的编号减1,循环到下次仍计算当前编号的向量(即原end)
                if(current != *end-1)
                {
                    //下面这个函数的方式应该是copy（from ,to， OK！
                    GCGE_ArrayCopyInSubspace(V+(*end -1)*ldV, V+current*ldV, nrows);
                }//end if (current != *end-1)
                (*end)--;
                current--;
                //如果这个时候需要输出提醒信息，就输出
                if(orth_para->print_orth_zero == 1)
                {
                    GCGE_Printf("In OrthonormalizationInSubspace, there is a zero vector!, "
                            "current = %d, start = %d, end = %d\n", current, start, *end);
                }//end if (orth_para->print_orth_zero == 1)
            }//end if (vout > orth_para->orth_zero_tol)
        }
    }
}//end for this subprogram

/** 
 *  V:         input,        要正交化的向量组
 *  start:     input,        要正交化的向量组在V中的起始位置
 *  end:       input|output, 要正交化的向量组在V中的终止位置
 *  B:         input,        正交化矩阵, 可以是NULL
 *  ops:       input,        操作
 *  para:      input,        参数
 *  workspace: input,        工作空间
 *
 *  块正交化过程：
 *  
 *  V = [V1, V2], V中有三个部分, 其中V1部分已经正交, 
 *  V1中有 size_V1(这里是nev) 个向量 
 *  V2中有 size_V2(这里是w_length=end-start) 个向量
 *
 *  临时存储空间使用workspace中的 evec(nev个), V_tmp(dim_xp-nev个), CG_p(1个)
 *
 *  下面的过程进行两次(相当于做一次重正交化)
 *  1. 去掉 V2 中 V1 方向的分量
 *     V2 = V2 - V1 * (V1^T * V2)
 *   > 计算 subspace_dtmp = V1^T * V2
 *          subspace_dtmp: GCGE_DOUBLE *数组
 *            size_V1 行, 即start行,w_length 列
 *   > 计算 V2 = V2 - V1 * subspace_dtmp
 *
 *  2. V2 自身做正交化 
 *     for current = start_V2:end_V2
 *         计算 d_tmp = (V[current:end_V2], V[current])
 *         norm = sqrt(d_tmp[0])
 *         如果norm > orth_zero_tol
 *            for V2_current = current:end_V2
 *                计算 V[V2_current] = V[V2_current] - d_tmp[V2_current-i]/norm * V[current]
 *            end
 *         否则，判定向量V[V2_current]为0
 *     end
 *
 */
void GCGE_SCBOrthonormalizationInSubspace(double *V, GCGE_INT ldV, GCGE_INT nrows, GCGE_INT start, 
        GCGE_INT *end, void *B, GCGE_INT ldB, GCGE_ORTH_PARA *orth_para, 
        GCGE_WORKSPACE *workspace, GCGE_OPS *ops)
{
    GCGE_Printf("Use GCGE_SCBOrthonormalizationInSubspace.\n");
    GCGE_INT    current = 0; //当前进行正交化操作的向量编号
    GCGE_INT    current_V2 = 0; //当前进行正交化操作的向量编号
    GCGE_INT    reorth_count = 0;
    GCGE_INT    orth_length = *end - start;
    GCGE_INT    old_orth_length = *end - start;
    GCGE_INT    orth_ncols = 1;
    GCGE_INT    max_reorth_count = orth_para->max_reorth_time;
    GCGE_DOUBLE orth_zero_tol = orth_para->orth_zero_tol;
    GCGE_DOUBLE alpha = 1.0;
    GCGE_DOUBLE beta = 0.0;
    GCGE_DOUBLE norm = 0.0;
    GCGE_DOUBLE test_temp = 0.0;
    GCGE_INT    test_i = 0;
    GCGE_DOUBLE *d_tmp = workspace->subspace_dtmp;
    GCGE_INT    tmp_eval_nonzero_start = 0;
    GCGE_INT    i = 0;

    if(start != 0)
    {
        for(reorth_count=0; reorth_count < max_reorth_count; reorth_count++)
        {
            orth_length = *end - start;
            //计算 d_tmp = V1^T * V2
            alpha = 1.0;
            beta  = 0.0;
            ops->DenseMatDotDenseMat("T", "N", &start, &orth_length,
                    &nrows, &alpha, V, &ldV, V+start*ldV, &ldV,
                    &beta, d_tmp, &start);
            //计算 V2 = V2 - V1 * subspace_dtmp
            alpha = -1.0;
            beta  = 1.0;
            ops->DenseMatDotDenseMat("N", "N", &ldV, &orth_length,
                    &start, &alpha, V, &ldV, d_tmp, &start,
                    &beta, V+start*ldV, &ldV);
        }
    }
    max_reorth_count = 1;
    for(reorth_count=0; reorth_count < max_reorth_count; reorth_count++)
    {
        //V2 自身做正交化 
        for(current = start; current < *end; ++current)
        {
            orth_length = *end - current;
            //计算 d_tmp = (V[current:end_V2], V[current])
            alpha = 1.0;
            beta  = 0.0;
            ops->DenseMatDotDenseMat("T", "N", &orth_length, &orth_ncols,
                    &nrows, &alpha, V+current*ldV, &ldV, V+current*ldV, &ldV,
                    &beta, d_tmp, &orth_length);
            //ops->DenseMatDotVec("T", &ldV, &orth_length, &alpha, 
            //        V+current*ldV, &ldV, V+current*ldV, &orth_ncols, 
            //        &beta, d_tmp, &orth_ncols);
            norm = sqrt(d_tmp[0]);
            if(norm > orth_zero_tol)
            {
                GCGE_ArrayScaleInSubspace(1.0/norm, V+current*ldV, nrows);
                for(current_V2=current+1; current_V2 < *end; current_V2++)
                {
                    if(fabs(d_tmp[current_V2 - current]/norm) > 10.0*orth_zero_tol)
                        GCGE_ArrayAXPBYInSubspace(-d_tmp[current_V2 - current]/norm, 
                            V+current*ldV, 1.0, V+current_V2*ldV, nrows);
                }
            }
            else
            {
                //如果v_current为0，就把最后一个向量拷贝到当前向量的位置
                //如果本向量就是最后一个向量，那就不拷贝(考虑如果是最后一个直接跳出)
                //同时(*end)--,把总向量的个数减1
                //同时current--,把当前向量的编号减1,循环到下次仍计算当前编号的向量(即原end)
                if(current != *end-1)
                {
                    //下面这个函数的方式应该是copy（from ,to， OK！
                    GCGE_ArrayCopyInSubspace(V+(*end -1)*ldV, V+current*ldV, nrows);
                }//end if (current != *end-1)
                (*end)--;
                current--;
                //如果这个时候需要输出提醒信息，就输出
                if(orth_para->print_orth_zero == 1)
                {
                    GCGE_Printf("In OrthonormalizationInSubspace, there is a zero vector!, "
                            "current = %d, start = %d, end = %d\n", current, start, *end);
                }//end if (orth_para->print_orth_zero == 1)
            }//end if (vout > orth_para->orth_zero_tol)
        }
    }
}//end for this subprogram

/*
 *  子空间正交化分批计算 GCGE_BlockOrthonormalizationInSubspace
 *  V:         input,        要正交化的向量组
 *  ldV:       input,        向量组V的leading dimension
 *  nrows:     input,        要进行正交化的向量长度
 *  end:       input|output, 要正交化的向量组在V中的终止位置
 *  orth_block_size: input,  每次正交化的向量个数
 *  ops:       input,        操作
 *  para:      input,        参数
 *  workspace: input,        工作空间
 *
 *  块正交化过程：
 *  
 *  V = [V1, V2, V3, ... , Vn], V中有n个部分, 
 *  n = (end+orth_block_size-1)/orth_block_size
 *  Vi中有 orth_size 个向量 
 *  orth_size = min(end-start, orth_block_size)
 *
 *  对 i = 0:n-1
 *  start = i * orth_block_size
 *  orth_size = min(end-start, orth_block_size)
 *  block_end = start + orth_size
 *
 *  1. 去掉 Vi 中 Vs = [V1, ... , V{i-1}] 方向的分量
 *     Vi = Vi - Vs * (Vs^T * Vi)
 *   > 计算 subspace_dtmp = Vs^T * Vi
 *          subspace_dtmp: GCGE_DOUBLE *数组
 *            size_Vs 行, 即 start 行, orth_size 列
 *   > 计算 Vi = Vi - Vs * subspace_dtmp
 *
 *  2. Vi 自身做正交化 
 *     for current = start:block_end
 *         计算 d_tmp = (V[current:block_end], V[current])
 *         norm = sqrt(d_tmp[0])
 *         如果norm > orth_zero_tol
 *            for Vi_current = current:block_end
 *                计算 V[Vi_current] = V[Vi_current] 
 *                - d_tmp[Vi_current-i]/norm * V[current]
 *            end
 *         否则，判定向量V[Vi_current]为0
 *     end
 *
 */
void GCGE_BlockOrthonormalizationInSubspace(GCGE_DOUBLE *V, GCGE_INT ldV, 
        GCGE_INT nrows, GCGE_INT *end, GCGE_INT orth_block_size,
        GCGE_OPS *ops, GCGE_PARA *para, GCGE_DOUBLE *subspace_dtmp)
{
    if(orth_block_size < 1)
    {
        orth_block_size = ((*end)/5 > 1) ? (*end)/5 : *end; 
    }  
    GCGE_INT num_block = ((*end)+orth_block_size-1)/orth_block_size;
    GCGE_INT current_block = 0;
    GCGE_INT start = 0;
    GCGE_INT orth_size = 0;
    GCGE_INT block_end = 0;
    GCGE_INT reorth_count = 0;
    GCGE_ORTH_PARA *orth_para = para->orth_para;
    GCGE_INT max_reorth_count = orth_para->max_reorth_time;

    GCGE_DOUBLE alpha = 0.0;
    GCGE_DOUBLE beta  = 0.0;

    for(current_block=0; current_block<num_block; current_block++)
    {
        //本批次要进行正交化的向量的起始列号,列数,与终点列号
        orth_size = ((*end)-start < orth_block_size)?((*end)-start):orth_block_size;
        block_end = start + orth_size;
        if(start > 0)
        {
            //减去前面的分量
            for(reorth_count=0; reorth_count < max_reorth_count; reorth_count++)
            {
                //计算 subspace_dtmp = Vs^T * Vi
                alpha = 1.0;
                beta  = 0.0;
                ops->DenseMatDotDenseMat("T", "N", &start, &orth_size,
                        &nrows, &alpha, V, &ldV, V+start*ldV, &ldV,
                        &beta, subspace_dtmp, &start);
                //计算 Vi = Vi - Vs * subspace_dtmp
                alpha = -1.0;
                beta  = 1.0;
                ops->DenseMatDotDenseMat("N", "N", &ldV, &orth_size,
                        &start, &alpha, V, &ldV, subspace_dtmp, &start,
                        &beta, V+start*ldV, &ldV);
            }
        }
        //自身正交化
        if(strcmp(para->x_orth_type, "bgs") == 0)
        {
            GCGE_BOrthonormalizationInSubspace(V+start*ldV, ldV, nrows, 0, &orth_size, 
                NULL, -1, orth_para, subspace_dtmp, ops);
        }
        else
        {
            GCGE_OrthonormalizationInSubspace(V+start*ldV, ldV, nrows, 0, &orth_size, 
                NULL, -1, orth_para);
        }
        //更新下一批次的start
        start += orth_size;
        if(start < block_end)
        {
            //说明自身正交化时有0向量出现
            //更新end
            *end -= block_end-start;
            //将最后几列向量拷贝到当前start位置
            memcpy(V+start*ldV, V+(*end)*ldV, (block_end-start)*ldV*sizeof(GCGE_DOUBLE));
        }
    }
}


//Multi Orthonormalization:
//先减去要正交化的向量中已正交化部分的分量, 
//再使用多重正交化方法进行自身的正交化
/*
 *  对向量组 V = [V1, V2] 进行正交化, 其中 V1 部分已经正交
 *  1. 去掉 V2 中 V1 方向的分量
 *     V2 = V2 - V1 * (V1^T * V2)
 *  2. V2 自身正交化
 *
 */
void GCGE_StableMultiOrthonormalization(void **V, GCGE_INT start, GCGE_INT *end, 
         void *B, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
    GCGE_INT mv_s[2];
    GCGE_INT mv_e[2];
    //GCGE_Printf("in GCGE_StableMultiOrthonormalization\n");

    GCGE_DOUBLE value_inner = 1.0;
    GCGE_DOUBLE norm_tmp = 0.0;
    GCGE_DOUBLE *subspace_dtmp = workspace->subspace_dtmp;
    void **V_tmp = workspace->V_tmp;
    GCGE_INT iter = 0;
    GCGE_INT i = 0;
    GCGE_DOUBLE Orth_Tol = para->orth_para->scbgs_reorth_tol;
    GCGE_INT max_reorth_time = para->orth_para->max_reorth_time;
    mv_s[0] = 0;
    mv_e[0] = start;
    mv_s[1] = start;
    mv_e[1] = *end;

    //先减去start之前的分量
    while(value_inner > Orth_Tol)
    {   
        iter ++;      

        GCGE_SubOrthonormalization(V, mv_s, mv_e, B, V_tmp, subspace_dtmp, ops);

        value_inner = 0.0;
        for(i=0;i<start*(*end - start);i++)
        { 
            norm_tmp = subspace_dtmp[i];
            value_inner += norm_tmp*norm_tmp;
        }
        value_inner = sqrt(value_inner);
        //GCGE_Printf("the %d-th orth, inner value=%e\n",iter, value_inner);    
        if(iter >= max_reorth_time)
            break;
    }//end while for value_inner

    //GCGE_Printf("in GCGE_StableMultiOrthonormalization, after GCGE_SubOrthonormalization\n");
    //然后start到end自身做正交化
    GCGE_MultiOrthonormalization(V, start, end, B, ops, para, workspace);

    //GCGE_Printf("line 1868, start: %d, end: %d\n",  start, *end);
}

//计算 V2 = V2 - V1 * (V1' * B * V2)
void *GCGE_SubOrthonormalization(void **V, GCGE_INT *start, GCGE_INT *end,
      void *B, void *V_tmp, GCGE_DOUBLE *subspace_dtmp, GCGE_OPS *ops)
{
    GCGE_INT length_1 = end[0]-start[0];
    GCGE_INT length_2 = end[1]-start[1];
    GCGE_INT mv_s[2];
    GCGE_INT mv_e[2];

    //计算 V_tmp = B * V2 
    if(B == NULL)
    {
        mv_s[0] = start[1];
        mv_e[0] = end[1];
        mv_s[1] = 0;
        mv_e[1] = length_2;
        ops->MultiVecAxpby(1.0, V, 0.0, V_tmp, mv_s, mv_e, ops);
    }
    else
    {
        mv_s[0] = start[1];
        mv_e[0] = end[1];
        mv_s[1] = 0;
        mv_e[1] = length_2;
        ops->MatDotMultiVec(B, V, V_tmp, mv_s, mv_e, ops);
    }

    //即计算 subspace_dtmp = V1^T * (B * V2)
    //计算 subspace_dtmp = V1^T * V_tmp
    mv_s[0] = start[0];
    mv_e[0] = end[0];
    mv_s[1] = 0;
    mv_e[1] = length_2;
    ops->MultiVecInnerProd(V, V_tmp, subspace_dtmp, "ns", mv_s, mv_e, length_1, ops);  

    //即计算 VV2 = V1 * (V1'*B*V2)
    //计算 V_tmp = V1 * subspace_dtmp
    mv_s[0] = start[0];
    mv_e[0] = end[0];
    mv_s[1] = 0;
    mv_e[1] = length_2;
    ops->MultiVecLinearComb(V, V_tmp, mv_s, mv_e, subspace_dtmp, length_1, NULL, -1, ops);
    
    //计算 V2 = V2 - V_tmp
    mv_s[0] = 0;
    mv_e[0] = length_2;
    mv_s[1] = start[1];
    mv_e[1] = end[1];
    ops->MultiVecAxpby(-1.0, V_tmp, 1.0, V, mv_s, mv_e, ops);
    
}

//对V(:, start:*end)自身进行正交化
/*
 *  使用递归的方式对 V 自身进行正交化
 *  将向量组 V 分成两部分, V = [V1, V2]
 *  1. 如果 V 中向量个数小于 max_direct_orth_length, 
 *     使用 BGS 直接对这部分进行正交化
 *     否则进行下面的过程
 *
 *  2. 调用 GCGE_MultiOrthonormalization 对 V1 部分进行正交化
 *  3. 减去 V2 中 V1 方向的分量
 *     V2 = V2 - V1 * (V1^T * V2)
 *  4. 调用 GCGE_MultiOrthonormalization 对 V2 部分进行正交化
 *
 */
void GCGE_MultiOrthonormalization(void **V, GCGE_INT start, GCGE_INT *end, void *B, 
      GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
    //GCGE_Printf("line 1873, start: %d, end: %d\n",  start, *end);
    GCGE_INT length = *end - start;
    //GCGE_Printf("in GCGE_MultiOrthonormalization, start: %d, end: %d, length: %d\n", 
    //          start, *end, length);
    //将start到end进行二分，并确定左右两侧的起始、终点位置，以及向量个数
    GCGE_INT length_1;
    GCGE_INT length_2;
    GCGE_INT mid;

    length_1 = length/2;
    length_2 = length - length_1;
    mid      = start + length_1;

    GCGE_INT       i = 0;
    GCGE_INT       mv_s[2];
    GCGE_INT       mv_e[2];
    GCGE_INT       iter = 0;
    GCGE_ORTH_PARA *orth_para = para->orth_para; 
    GCGE_DOUBLE    Orth_Tol = orth_para->scbgs_reorth_tol;
    GCGE_DOUBLE    *subspace_dtmp = workspace->subspace_dtmp;
    GCGE_DOUBLE    value_inner = 1.0;
    GCGE_DOUBLE    norm_tmp = 0.0;
    void           **V_tmp = workspace->V_tmp;
    void           *vec;
    void           *vec_tmp;

    GCGE_INT max_direct_orth_length = para->orth_para->max_direct_orth_length;
    if(length <= max_direct_orth_length)
    {
        //使用GCGE_SubOrthonormalizationSelfWithEigen这种正交化方式，
        //精度要求较高时，很容易出问题
        GCGE_SubOrthonormalizationSelfBGS(V, start, end, B, 
              para, ops, workspace);
    }
    else
    {

        //左侧
        //因为可以直接处理length==1, 所以直接调用GCGE_MultiOrthonormalization
        GCGE_MultiOrthonormalization(V, start, &mid, B, ops, para, workspace);
        
        //右侧
        //如果左侧有0向量出现，先将右侧最后的向量拷贝到前面
        if(mid + length_2 < *end)
        {
            GCGE_INT new_end_2 = mid + length_2;
            GCGE_INT n_zero = *end - new_end_2;
            mv_s[0] = mid;
            mv_e[0] = mid + n_zero;
            mv_s[1] = new_end_2;
            mv_e[1] = *end;
            ops->MultiVecSwap(V, V, mv_s, mv_e, ops);
            *end = new_end_2;
        }
	if(mid-start > 0)
	{
            //右侧减去左侧的分量
            while(value_inner > Orth_Tol)
            {   
                iter ++;      
                mv_s[0] = start;
                mv_e[0] = mid;
                mv_s[1] = mid;
                mv_e[1] = *end;
                GCGE_SubOrthonormalization(V, mv_s, mv_e, B, V_tmp, subspace_dtmp, ops);
                value_inner = 0.0;
                for(i=0;i<length_1*length_2;i++)
                { 
                    norm_tmp = subspace_dtmp[i];
                    value_inner += norm_tmp*norm_tmp;
                }
                value_inner = sqrt(value_inner);
                //GCGE_Printf("the %d-th orth, inner value=%e\n",iter, value_inner);    
                if(iter >= orth_para->max_reorth_time)
                    break;
            }//end while for value_inner
	    if(value_inner > Orth_Tol)
	    {
	       GCGE_Printf("in GCGE_SubOrthonormalization, V1: %d-%d, V2: %d-%d, reorth_count: %d, value_inner: %e\n", 
		     start, mid-1, mid, *end, orth_para->max_reorth_time, value_inner);
	    }
	}
        
        GCGE_MultiOrthonormalization(V, mid, end, B, ops, para, workspace);
    }
}

void GCGE_SubOrthonormalizationSelfBGS(void **V, GCGE_INT start, GCGE_INT *end, 
      void *B, GCGE_PARA *para, GCGE_OPS *ops, GCGE_WORKSPACE *workspace)
{
    GCGE_INT current;
    void **V_tmp = workspace->V_tmp;
    void *vec_current;
    void *vec_tmp;
    GCGE_INT j = 0;
    GCGE_INT mv_s[2];
    GCGE_INT mv_e[2];
    GCGE_INT length = 0;
    GCGE_DOUBLE *d_tmp = workspace->subspace_dtmp;
    GCGE_DOUBLE norm_value = 0.0;
    GCGE_DOUBLE orth_zero_tol = para->orth_para->orth_zero_tol;
    //现在对W本身进行正交化 
    for(current=start;current<*(end);current++)
    {
        //vec_current = W(:,current)
        ops->GetVecFromMultiVec(V, current, &vec_current);
        ops->GetVecFromMultiVec(V_tmp, 0, &vec_tmp);
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
        mv_s[0] = current;
        mv_e[0] = *end;
        mv_s[1] = 0;
        mv_e[1] = 1;
        //d_tmp = vec_tmp*V(:,start:end)
        length = mv_e[0] - mv_s[0];
        //统一做内积
        ops->MultiVecInnerProd(V, V_tmp, d_tmp, "ns", mv_s, mv_e, length, ops);    
        //V[current]的向量范数
        norm_value = sqrt(d_tmp[0]);      
        //下面看看目前处理的向量是否是一个零向量          
        if(norm_value > orth_zero_tol)
        {
            //如果vec_current不为0，就做归一化
            ops->GetVecFromMultiVec(V, current, &vec_current);
            //vec_current = 1/norm_value * vec_current
            ops->VecAxpby(0.0, vec_current, 1.0/norm_value, vec_current);    
            ops->RestoreVecForMultiVec(V, current, &vec_current);
            for(j=current+1;j<*end;j++)
            {
                ops->GetVecFromMultiVec(V, current, &vec_current);
                ops->GetVecFromMultiVec(V, j, &vec_tmp);
                //vec_tmp = vec_tmp - d_tmp[j-current]/norm_value *vec_current
                ops->VecAxpby(-d_tmp[j-current]/norm_value, vec_current, 1.0, vec_tmp);
                ops->RestoreVecForMultiVec(V, j, &vec_tmp);
                ops->RestoreVecForMultiVec(V, current, &vec_current);
            }//end for(j=current+1;j<*(end);j++)
        }
        else
        {
            //如果vec_current为0，就把此向量与最后一个向量进行指针交换
            //同时(*end)--,把总向量的个数减1
            //同时current--,把当前向量的编号减1,循环到下次仍计算当前编号的向量(即原end)
            mv_s[0] = current; //表示swap中第一个向量组的起始位置
            mv_s[1] = (*end)-1; //表示swap中第二个向量组的起始位置
            mv_e[0] = current + 1; //表示swap中第一个向量组的终止位置
            mv_e[1] = *end; //表示swap中第二个向量组的终止位置
            //向量移动（由用户提供），目的: V(:,current:end-1) = V(:, )
            ops->MultiVecSwap(V, V, mv_s, mv_e, ops);
            (*end)--;
            current--;
            if(para->orth_para->print_orth_zero == 1)
            {
                GCGE_Printf("In GCGE_SubOrthonormalizationSelfBGS, there is a zero vector!, "
                        "current = %d, start = %d, end = %d\n", current, start, *end);              
            }//end if(orth_para->print_orth_zero == 1)        
        }//end if(vout > orth_para->orth_zero_tol)      
    }//end for(current=dim_xp;current<*(end);current++)      
}

/*
 *  使用递归的方式对 V 自身进行正交化
 *  将向量组 V 分成两部分, V = [V1, V2]
 *  1. 如果 V 中向量个数小于 max_direct_orth_length, 
 *     使用 BGS 直接对这部分进行正交化
 *     否则进行下面的过程
 *
 *  2. 调用 GCGE_MultiOrthonormalizationInSubspace 对 V1 部分进行正交化
 *  3. 减去 V2 中 V1 方向的分量
 *     V2 = V2 - V1 * (V1^T * V2)
 *  4. 调用 GCGE_MultiOrthonormalizationInSubspace 对 V2 部分进行正交化
 *
 */
void GCGE_MultiOrthonormalizationInSubspace(double *V, GCGE_INT ldV, GCGE_INT nrows, 
        GCGE_INT start, GCGE_INT *end, void *B, GCGE_INT ldB, 
        GCGE_ORTH_PARA *orth_para, GCGE_WORKSPACE *workspace, GCGE_OPS *ops)
{
    //将start到end进行二分，并确定左右两侧的起始、终点位置，以及向量个数
    GCGE_INT length   = *end - start;
    GCGE_INT length_1 = length/2;
    GCGE_INT length_2 = length - length_1;
    GCGE_INT mid      = start + length_1;

    GCGE_INT    reorth_count = 0;
    GCGE_INT    max_reorth_count = orth_para->max_reorth_time;
    GCGE_DOUBLE alpha = 1.0;
    GCGE_DOUBLE beta = 0.0;
    GCGE_DOUBLE *subspace_dtmp = workspace->subspace_dtmp;
    GCGE_INT    max_direct_orth_length = orth_para->max_direct_orth_length;
    GCGE_INT    old_length = length;
    GCGE_DOUBLE value_tmp;

    //GCGE_Printf("line 2389, start: %d, end: %d, length: %d\n", 
    //              start, *end, length);
    //if(length <= max_direct_orth_length)
    if(length == 1)
    {
        //GCGE_OrthonormalizationInSubspace(V+start*ldV, ldV, nrows, 0, &length, 
        //      B, ldB, orth_para);
        value_tmp = GCGE_ArrayNormInSubspace(V+start*ldV, nrows);
        if(value_tmp > orth_para->orth_zero_tol)
        {
            GCGE_ArrayScaleInSubspace(1.0/value_tmp, V+start*ldV, nrows);
        }//做完归一化的情况，下面考虑是 0 向量的情况
        else
        {
            *end = start;
        }
        //value_tmp = GCGE_ArrayNormInSubspace(V+start*ldV, nrows);
        //GCGE_Printf("line 2406, start: %d, end: %d, length: %d, ip: %e\n", 
        //      start, *end, length, value_tmp);
        //if(length < old_length)
        //    *end = start + length;
    }
    else
    {
        GCGE_MultiOrthonormalizationInSubspace(V, ldV, nrows, start, &mid, B, ldB, 
              orth_para, workspace, ops);
        if(mid + length_2 < *end)
        {
            GCGE_INT new_end_2 = mid + length_2;
            GCGE_INT n_zero = *end - new_end_2;
            length_1 -= n_zero;
            memcpy(V+mid*ldV, V+new_end_2*ldV, n_zero*ldV*sizeof(GCGE_DOUBLE));
            *end = new_end_2;
        }
	if(length_1 > 0)
	{
            for(reorth_count=0; reorth_count < max_reorth_count; reorth_count++)
            {
                //计算 d_tmp = V1^T * V2
                alpha = 1.0;
                beta  = 0.0;
                ops->DenseMatDotDenseMat("T", "N", &length_1, &length_2,
                        &nrows, &alpha, V+start*ldV, &ldV, V+mid*ldV, &ldV,
                        &beta, subspace_dtmp, &length_1);
                //计算 V2 = V2 - V1 * subspace_dtmp
                alpha = -1.0;
                beta  = 1.0;
                ops->DenseMatDotDenseMat("N", "N", &ldV, &length_2, &length_1, 
                        &alpha, V+start*ldV, &ldV, subspace_dtmp, &length_1,
                        &beta, V+mid*ldV, &ldV);
            }
	}
        GCGE_MultiOrthonormalizationInSubspace(V, ldV, nrows, mid, end, B, ldB, 
              orth_para, workspace, ops);
    }
}

//先减去要正交化的向量中已正交化部分的分量, 
//再使用多重正交化方法进行自身的正交化
/*
 *  对向量组 V = [V1, V2] 进行正交化, 其中 V1 部分已经正交
 *  1. 去掉 V2 中 V1 方向的分量
 *     V2 = V2 - V1 * (V1^T * V2)
 *  2. V2 自身正交化
 *
 */
void GCGE_StableMultiOrthonormalizationInSubspace(double *V, GCGE_INT ldV, 
       GCGE_INT nrows, GCGE_INT start, GCGE_INT *end, void *B, GCGE_INT ldB, 
       GCGE_ORTH_PARA *orth_para, GCGE_WORKSPACE *workspace, GCGE_OPS *ops)
{
    GCGE_INT    reorth_count = 0;
    GCGE_INT    orth_length = *end - start;
    GCGE_INT    old_orth_length = *end - start;
    GCGE_INT    max_reorth_count = orth_para->max_reorth_time;
    GCGE_DOUBLE alpha = 1.0;
    GCGE_DOUBLE beta = 0.0;
    GCGE_DOUBLE *d_tmp = workspace->subspace_dtmp;
    if(start != 0)
    {
        for(reorth_count=0; reorth_count < max_reorth_count; reorth_count++)
        {
            orth_length = *end - start;
            //计算 d_tmp = V1^T * V2
            alpha = 1.0;
            beta  = 0.0;
            ops->DenseMatDotDenseMat("T", "N", &start, &orth_length,
                    &nrows, &alpha, V, &ldV, V+start*ldV, &ldV,
                    &beta, d_tmp, &start);
            //计算 V2 = V2 - V1 * subspace_dtmp
            alpha = -1.0;
            beta  = 1.0;
            ops->DenseMatDotDenseMat("N", "N", &ldV, &orth_length,
                    &start, &alpha, V, &ldV, d_tmp, &start,
                    &beta, V+start*ldV, &ldV);
        }
    }
    GCGE_MultiOrthonormalizationInSubspace(V, ldV, nrows, start, end, B, ldB, 
        orth_para, workspace, ops);
}

