/*
 * =====================================================================================
 *
 *       Filename:  gcge_eigsol.c
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

#include "gcge_eigsol.h"



/**
 * @brief 计算 A * evec = eval * B * evec 
 *    
 *    Initialization: 
 *       对evec设定随机初值, 把evec复制给X
 *       对X正交化,在X张成的子空间中计算RR问题
 *       
 *    Loop: 
 *       对未收敛的X计算P
 *       对未收敛的X计算W
 *       对W进行正交化
 *       在[X P W]张成的子空间中计算RR问题
 *       检查特征值重数，调整X的长度
 *       检查收敛性, 得到block_size, unlock, num_unlock
 *
 *    Finalization:
 *       把X复制给evec
 *       把subspace_eval复制给eval
 *
 * @param A
 * @param B
 * @param eval
 * @param evec
 * @param nev  要求的特征值个数
 * @param ops
 */
void GCGE_Default_SolveEigenvalues(void *A, void *B, GCGE_DOUBLE *eval, void **evec, GCGE_PARA *para, GCGE_OPS *ops)
{
	GCGE_WORKSPACE *workspace;
	GCGE_WORKSPACE_Create(&workspace, para, ops, A);

	GCGE_INT    max_size_X       = workspace->space_size_RitzVec;      /* 最大的size_X, 即最多求解Ritz值的个数  */
	GCGE_INT    max_iter_num     = para->ev_max_iter;    /* 最大迭代次数 */

	GCGE_INT    *unlock          = workspace->unlock;          /* 存储未收敛的特征对编号, 分配的长度为? */
	GCGE_DOUBLE *subspace_matrix = workspace->subspace_matrix; /* V^T A V,                分配的长度为? */
	GCGE_DOUBLE *subspace_eval   = workspace->subspace_eval;   /* 子空间特征向量,         分配的长度为? */
	GCGE_DOUBLE *subspace_evec   = workspace->subspace_evec;   /* 子空间特征向量,         分配的长度为? */
	void        **V              = workspace->V;               /* 存储 [X P W],           分配的长度为? */
	void        **RitzVec        = workspace->RitzVec;         /* 存储Ritz向量,           分配的长度为? */
    void        *Orth_mat; //Orth_mat正交化时计算内积所用的矩阵: B 或者 A，通常是B
    void        *RayleighRitz_mat; //RayleighRitz_mat计算子空间矩阵时所用矩阵: A或者B,通常用A

	GCGE_INT    size_X , size_P , size_W ; /* X P W 的长度, start_X == 0 */
	GCGE_INT    start_X, start_P, start_W; /* X P W 在V中的起始列号      */
	GCGE_INT    end_X  , end_P  , end_W  ; /* X P W 在V中的结束列号+1    */
	GCGE_INT    size_V_old, size_X_old;

	GCGE_INT    block_size = 0;            /* 分批计算特征对的个数       */
	GCGE_INT    num_iter   = 0;            /* 迭代次数                   */
	GCGE_INT    mv_s[2]    = {0, 0};
	GCGE_INT    mv_e[2]    = {0, 0};

    //如果使用B内积，那么用矩阵B进行正交化，用矩阵A计算子空间矩阵
    //如果使用A内积，那么用矩阵A进行正交化，用矩阵B计算子空间矩阵
    //strcmp函数：当两个字符串相同时返回0, 如果orth_type==B, 返回0,运行else
    if(strcmp(para->orth_type, "B"))
    {
        //这个是表示不是用矩阵B做内积的情况， 因为 strcmp(para->orth_type, "B")不为零的时候表示的
        //orth_type不是B
        Orth_mat = A;
        RayleighRitz_mat = B;
    }
    else
    {   //表示的使用矩阵B做内积的情况
        Orth_mat = B;
        RayleighRitz_mat = A;
    }

	/* Initialization start 
	 *
	 * 对X做正交化
	 */
	size_X  = nev;  size_P = 0;  size_W = 0;
	start_X = 0;     end_X = size_X;
	start_P = end_X; end_P = start_P+size_P;
	start_W = end_P; end_W = start_W+size_W;
	size_V  = end_W;
	/*
	 * 如果用户不给初值，这里给随机初值，用户需要提供给向量设随机值的函数
	 */
	if(workspace->given_init_evec == 0)
	{
		mv_s[0] = start_X; mv_e[0] = end_X;
		ops->MultiVecSetRandomValue(V, mv_s, mv_e, ops);
	}
	else
	{
		/* 把用户提供的evec初值copy给V */
		mv_s[0] = 0;       mv_e[0] = nev; 
		mv_s[1] = start_X; mv_e[1] = end_X;
		ops->MultiVecAxpby(1.0, evec, 0.0, V, mv_s, mv_e, ops);
	}
	/* 
	 * 对初始近似特征向量做B正交化, 会删去线性相关的向量
	 * 返回end_X
	 */

    GCGE_Orthogonal(V, start_X, &end_X, Orth_mat, ops, para->orth_para, 
			workspace->V_tmp, workspace->subspace_dtmp);

	size_X = end_X - start_X;
	if(size_X < nev)
	{
		printf("The initial eigenvectors is linearly dependent!\n");
	}
	/* 子空间矩阵的大小 */
	size_V = end_X;

	/* Initialization end */

	/* LOOP start 
	 *
	 * size_V是上一次迭代[X P W]的长度, end_X是上一次的X长度
	 * 基于上次迭代的[X P W]更新P, 再利用RitzVec更新X, 利用已更新的X更新W 
	 * 利用新的[X P W] 更新subspace_matrix, 并求特征值和特征向量subspace_eval, subspace_evec
	 * 检查重数更新size_X, 并得到RitzVec
	 * 再检查收敛性更新unlock
	 */
	num_unlock = nev;
	num_iter   = 0;
	while((num_unlock > 0)&&(num_iter < max_iter_num))
	{
		/*
		 * 若size_V == end_X, 表明XP和XW都没有，那么不能计算P, 直接返回start_P = end_X; end_P = start_P;
		 */
		if (size_V == end_X)
		{
			start_P = end_X;
			end_P   = start_P;
		}
		else 
		{
			/*
			 * size_X 会被上一次迭代最后检查重数而更新, 它是Ritz向量的个数, 
			 * 此时end_X只是上一次的X个数, 即当前V中X的长度, 
			 *
			 * 计算P=X_{j}-X_{j-1}, 是通过subspace的正交操作实现P的计算以及对X的正交化
			 * subspace_evec
			 *    XX  PX                 XX  O
			 *    XP  PP  => 复制并置零  XP  XP  => 正交化
			 *    XW  PW                 XW  XW
			 *
			 * XX的行列数分别是end_X和size_X, 
			 * PX的行列数分别是end_X和block_size, 
			 *
			 * 对subspace_evec, 将unlock[j]列进行复制，并对前end_X行置零
			 * 进行subspace上的正交化之后
			 * 返回start_P = size_X, end_P = start_P + block_size ... (若没有线性相关出现)
			 * 并做线性组合, 将P更新到V的从start_P到end_P-1列
			 *
			 * 返回start_P, end_P
			 */
			GCGE_ComputeP(V, start_P, &end_P, subspace_evec, end_V, unlock, 
				size_X, size_X_old, workspace->V_tmp, para, ops);
		}
		size_P = end_P - start_P;

		/* 当第一次迭代时，只计算X生成的subspace的特征值问题 */
		start_W = end_P; 
		if (num_iter == 0)
		{
			end_W = start_W;
		}
		else
		{
			/*
			 * 将V中的X更新, 即V中的X部分与RitzVec交换指针
			 * 这里之所以再次赋值end_X是因为endX是上次迭代时V中的X长度
			 */
			end_X   = firstX + size_X;
			mv_s[0] = fisrtX; mv_e[0] = end_X;
			mv_s[1] = 0;      mv_e[1] = size_X;
			ops->MultiVecSwap(V, RitzVec, mv_s, mv_e, ops);
			/*
			 * 计算AW = \lambda BX, 只计算第unlock[j]列的X对应的W, 0<=j<block_size
			 * 并将W对[X, P]做正交化, 
			 * 返回end_W, 即size_W = end_W-start_W. 这是由于过程中会将线性相关的部分去掉
			 * 若没有出现线性相关，那么size_W == block_size
			 */
			size_W = unconv_bs;
			end_W  = start_W + size_W;
			GCGE_ComputeW(A, B, V, start_W, size_W, subspace_eval, ops, para, workspace);
			GCGE_Orthogonal(V, start_W, &end_W, Orth_mat, ops, para->orth_para, 
				workspace->V_tmp, workspace->subspace_dtmp);
		}
		size_W = end_W - start_W;
		/*
		 * 计算子空间矩阵
		 * subspace_matrix = [X, P, W]^T*A*[X, P, W]
		 * subspace_matrix是size_V阶方阵，且leading dimension也是size_V
		 */
		size_V_old = size_V;
		size_V     = end_W;

		ComputeSubspaceMatrix(A, V, subspace_matrix, size_V, ops);
		GCGE_ComputeSubspaceMatrix(subspace_matrix, size_V, ldm_old, RayleighRitz_mat,
				V, size_V, subspace_evec, end_P, ops, 
				workspace->subspace_dtmp, workspace->V_tmp);
		/*
		 * 求解size_X个计算子空间矩阵的特征对
		 * 存储在subspace_eval, subspace_evec
		 */
		GCGE_ComputeSubspaceEigenpairs(subspace_matrix, size_V, subspace_eval, 
				subspace_evec, ops, para, workspace);
		/*
		 * 检查特征值重数, 确定size_X, 即本次更新X的个数, 
		 * 完成检查重数之后, 会更改计算Ritz向量的个数, 
		 * 但end_X不更新, 存储的是上一次的size_X 
		 */
		size_X_old = size_X;
		GCGE_CheckEvalMultiplicity(nev+2, max_size_X, &size_X, subspace_eval);
		/*
		 * 用子空间特征向量与基向量V计算size_X个Ritz向量, 即
		 * RitzVec[j] = [X P W] subspace_evec[j], 0<=j<size_X
		 */
		GCGE_ComputeX(V, size_V, subspace_evec, RitzVec, size_X, ops);
		/*
		 * 检查收敛性, 计算特征对的残差，更新block_size, num_unlock和unlock
		 *
		 * nev是要算的特征值个数，size_X是计算出的特征值个数
		 * TODO 需要加入形参tol用以检查收敛性
		 * 
		 * 计算
		 * idx = unlock[j]  0<=j<block_size
		 * A RitzVec[ idx ] - subspace_eval[ idx ] B RitzVec[ idx ] 的模并判断是否收敛
		 *
		 * 返回block_size, 即下次迭代时, 对unlock[j], 0<=j<block_size的X计算P和W
		 * 返回unlock和num_unlock, 记录所有未收敛的特征对编号
		 */
		GCGE_CheckConvergence(nev, A, B, 
				size_X, subspace_eval, RitzVec, 
				&block_size, &num_unlock, unlock, ops);
		for(i=0; i<nev; i++)
		{
			printf("num_iter: %d, num_unlock: %d, eval[%d]: %e\n", 
					num_iter, num_unlock, i, subspace_eval[i]);
		}
		++num_iter;
	}
	/* LOOP end */

	/* Finalization start 
	 *
	 * 把计算得到的近似特征对拷贝给eval, evec
	 *
	 * TODO 这里可以用BLAS进行copy特征值
	 *      并且要检查是否已经收敛了nev个特征值
	 *      考虑nev是否应该是一个返回值
	 */
	memcpy(eval, subspace_eval, nev*sizeof(GCGE_DOUBLE));
	mv_s[0] = 0; mv_e[0] = nev;
	mv_s[1] = 0; mv_e[1] = nev;
	ops->MultiVecAxpby(1.0, RitzVec, 0.0, evec, mv_s, mv_e, ops);
	GCGE_WORKSPACE_Free(&workspace, ops);
	/* Finalization end */
}




/*  brief: 检查(eval, X)的收敛性,返回收敛信息num_unlock,unconv_bs,unlock
 *     
 *     对unlock中的num_unlock个编号的特征向量计算残差，判断收敛性
 *     允许收敛的特征对不连续
 *     一旦出现连续10个特征对未收敛，就默认后面的都未收敛，不全部检测
 *      
 *     需要的参数：unlock (in \ out) : 未收敛的编号
 *                 num_unlock (in \ out) : 未收敛的特征对个数
 *                 unconv_bs (out) : 本批次中未收敛的特征对个数
 *                 res (out) : 记录残差
 *
 */
//计算残差，检查收敛性，并获取未收敛的特征对编号及个数
//A,B是用来计算残差的矩阵(GCGE_MATRIX), evec是特征向量(GCGE_Vec)
//unlock中所包含的未收敛的特征值有可能是为了保证收敛速度而额外算得，所以要跟nev的关系弄清楚
void GCGE_CheckConvergence(void *A, void *B, GCGE_INT *num_unlock, GCGE_DOUBLE *eval, void **evec, 
        GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace)
{
    /* TODO */
    GCGE_INT    i, j, k, *unlock = workspace->unlock,
                unconv_bs = 0, max_conv_idx = 0,
                max_ind = 0, min_ind = 0;
    GCGE_DOUBLE res_norm, evec_norm, residual, 
                max_res = 0.0, min_res = 0.0, sum_res = 0.0,
                ev_tol = para->ev_tol,
                *res = para->res;
    char        *conv_type = para->conv_type;
    void        *tmp1, *tmp2, *ev;
    ops->GetVecFromMultiVec(workspace->V_tmp, 0, &tmp1);
    ops->GetVecFromMultiVec(workspace->V_tmp, 1, &tmp2);
    for( j=0; j<(*num_unlock); j++ )
    {
        i = unlock[j];
        ops->GetVecFromMultiVec(evec, i, &ev);
        //计算残量
        ops->MatDotVec(A, ev, tmp1);
        if(B == NULL)
        {
            ops->VecAxpby(-eval[i], ev, 1.0, tmp1);
        }
        else
        {
            ops->MatDotVec(B, ev, tmp2);
            ops->VecAxpby(-eval[i], tmp2, 1.0, tmp1);
        }
        //tmp1是残差向量
        //GCGE_VecNorm：这个函数可以放到ops中去。
        res_norm  = GCGE_VecNorm(tmp1, ops);
        //printf("eigenvalue = %5.15f,  i = %d,  res_norm = %5.15f\n", eval[i], i, res_norm);
        evec_norm = GCGE_VecNorm(ev, ops);
        //printf("||evec[:,i]||=%5.15f\n",evec_norm);
          //计算相对/绝对残差
        if(strcmp(conv_type, "R"))
        {
          residual  = res_norm/evec_norm;
        }
        else
        {
             //计算相对残差的时候需要再仔细考虑一下
           residual  = res_norm/evec_norm/(1.0>fabs(eval[i])?1.0:fabs(eval[i]));
        }
        res[i] = residual;
        //统计最大、最小、总残差
        sum_res += residual;
        if(j == 0)
        {
            max_res = residual;
            min_res = residual;
        }
        else
        {
            if(residual > max_res)
            {
                max_res = residual;
                max_ind = i;
            }
            if(residual < min_res)
            {
                min_res = residual;
                min_ind = i;
            }
        }//end if(j==0)
        //统计未收敛的特征对编号unlock及数量
        if(residual > ev_tol)
          {
            unlock[unconv_bs] = i;
            unconv_bs += 1;
            //只比收敛的最后一个多算8个
            //在unlock数组中检查到上一次收敛位置往后8位
            //这里的 8 最好是在para中设置一下
          if(j > max_conv_idx+8)
            {
                for(k=j+1; k<(*num_unlock); k++)
                {
                    unlock[unconv_bs] = unlock[k];
                    unconv_bs += 1;
                }
             j = (*num_unlock)+1; //直接退出 for j的循环
          }//end  if(j > max_conv_idx+8)            
          }
        else
        {
            max_conv_idx = j;
        }
        ops->RestoreVecForMultiVec(evec, i, &ev);
    }//end for j
    ops->RestoreVecForMultiVec(workspace->V_tmp, 0, &tmp1);
    ops->RestoreVecForMultiVec(workspace->V_tmp, 1, &tmp2);
     //整个0到nev没有收敛特征对的个数
    *num_unlock = unconv_bs;
     //当前的BLOCK中没有收敛的个数
    workspace->unconv_bs = (unconv_bs < para->block_size) ? unconv_bs : para->block_size;

    para->max_res = max_res;
    para->max_ind = max_ind;
    para->min_res = min_res;
    para->min_ind = min_ind;
    para->sum_res = sum_res;
    //printf("check : sum_res = %e\n",para->sum_res);

    para->num_iter += 1;
    //打印收敛信息
    /*
    if(para->print_result == 1)
    {
#if USE_MPI
        GCGE_PrintConvInfoParallel(solver);
#else
        GCGE_PrintConvInfo(solver);
#endif
    }
    */

#if 0 
    GCGE_DOUBLE t2 = GCGE_GetTime();
    para->stat_para->ite_time = t2 - para->stat_para->ite_time;
    para->stat_para->PartTimeTotal->Conv_Time += t2-t1;
    para->stat_para->PartTimeInOneIte->Conv_Time = t2-t1;
#endif
}//end for this subprogram

/* brief:  从第nev个特征值开始检查重数，如果后一个特征值与前一个的比值在0.8到1.2之间，
 *         就认为可能是重特征值
 *
 *  需要的参数：需要检查的特征值的起始和终点位置
 *              start (in)
 *              end   (in \ out)
 */
//检查特征值重数，更新dim_x
void GCGE_CheckEvalMultiplicity(GCGE_INT start,GCGE_INT nev, GCGE_INT end, 
                                GCGE_INT *dim_x, GCGE_DOUBLE *eval)
{
    GCGE_INT tmp, i;
    tmp = start;
     if(start < nev)
     {
        printf("The parameter setting is wrong!\n");
        exit(0);
      }  
     
    //dsygv求出的特征值已经排序是ascending,从小到大
    //检查特征值的数值确定下次要进行计算的特征值个数
    for( i=start; i<end; i++ )
    {
        if((fabs(fabs(eval[tmp]/eval[nev])-1))<0.2)
            tmp += 1;
        else
            break;
    }
    //更新新的X的向量个数
    *dim_x = tmp;
}

