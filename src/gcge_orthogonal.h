/*
 * =====================================================================================
 *
 *       Filename:  gcge_orthogonal.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2018年09月24日 09时50分16秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#ifndef _GCGE_ORTHOGONAL_H_
#define _GCGE_ORTHOGONAL_H_

#include <string.h>
#include "gcge_config.h"

#include "gcge_para.h"
#include "gcge_ops.h"

#include "gcge_workspace.h"
#include "gcge_utilities.h"

//GCGE算法中用到的正交化操作
/**
 * @brief 
 *
 * @param V
 * @param start
 * @param[in|out] end
 * @param B
 * @param ops
 * @param orth_para
 * @param workspace
 */
void GCGE_Orthonormalization(void **V, GCGE_INT start, GCGE_INT *end, 
      void *B, GCGE_OPS *ops, GCGE_PARA *para, 
      void **V_tmp, GCGE_DOUBLE *d_tmp);

void GCGE_OrthonormalizationInSubspace(double *V, GCGE_INT ldV, GCGE_INT nrows, GCGE_INT start, GCGE_INT *end, 
      void *B, GCGE_INT ldB, GCGE_ORTH_PARA *orth_para);
void GCGE_BOrthonormalization(void **V, GCGE_INT start, GCGE_INT *end, 
      void *B, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace);
void GCGE_CBOrthonormalization(void **V, GCGE_INT start, GCGE_INT *end, 
      void *B, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace);
void GCGE_SCBOrthonormalization(void **V, GCGE_INT start, GCGE_INT *end, 
      void *B, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace);

void GCGE_BOrthonormalizationInSubspace(double *V, GCGE_INT ldV, GCGE_INT nrows, GCGE_INT start, 
      GCGE_INT *end, void *B, GCGE_INT ldB, GCGE_ORTH_PARA *orth_para, 
      GCGE_DOUBLE *d_tmp, GCGE_OPS *ops);
void GCGE_SCBOrthonormalizationInSubspace(double *V, GCGE_INT ldV, GCGE_INT nrows, GCGE_INT start, 
      GCGE_INT *end, void *B, GCGE_INT ldB, GCGE_ORTH_PARA *orth_para, 
      GCGE_WORKSPACE *workspace, GCGE_OPS *ops);
void GCGE_BlockOrthonormalizationInSubspace(GCGE_DOUBLE *V, GCGE_INT ldV, 
        GCGE_INT nrows, GCGE_INT *end, GCGE_INT orth_block_size,
        GCGE_OPS *ops, GCGE_PARA *para, GCGE_DOUBLE *subspace_dtmp);
void GCGE_SCBOrth_Self(void **V, GCGE_INT start, GCGE_INT *end, 
              void *B, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace);
void GCGE_SCBOrth_Minus(void **V, GCGE_INT start, GCGE_INT *end, 
        void *B, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace);

//Multi Orthonormalization:
void GCGE_StableMultiOrthonormalization(void **V, GCGE_INT start, GCGE_INT *end, 
         void *B, GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace);
void *GCGE_SubOrthonormalization(void **V, GCGE_INT *start, GCGE_INT *end,
      void *B, void *V_tmp, GCGE_DOUBLE *subspace_dtmp, GCGE_OPS *ops);
void GCGE_MultiOrthonormalization(void **V, GCGE_INT start, GCGE_INT *end, void *B, 
      GCGE_OPS *ops, GCGE_PARA *para, GCGE_WORKSPACE *workspace);
void GCGE_SubOrthonormalizationSelfBGS(void **V, GCGE_INT start, GCGE_INT *end, 
      void *B, GCGE_PARA *para, GCGE_OPS *ops, GCGE_WORKSPACE *workspace);
void GCGE_MultiOrthonormalizationInSubspace(double *V, GCGE_INT ldV, GCGE_INT nrows, 
        GCGE_INT start, GCGE_INT *end, void *B, GCGE_INT ldB, 
        GCGE_ORTH_PARA *orth_para, GCGE_WORKSPACE *workspace, GCGE_OPS *ops);
void GCGE_StableMultiOrthonormalizationInSubspace(double *V, GCGE_INT ldV, 
       GCGE_INT nrows, GCGE_INT start, GCGE_INT *end, void *B, GCGE_INT ldB, 
       GCGE_ORTH_PARA *orth_para, GCGE_WORKSPACE *workspace, GCGE_OPS *ops);
#endif
