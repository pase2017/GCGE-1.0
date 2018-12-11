/*
 * =====================================================================================
 *
 *       Filename:  gcge_utilities.h
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
 /*  */
#ifndef  _GCGE_UTILITIES_H_
#define  _GCGE_UTILITIES_H_

#include "gcge_config.h"
#include <stdarg.h>
#if GCGE_USE_MPI
#include <mpi.h>
#endif

GCGE_DOUBLE GCGE_GetTime();
GCGE_INT    GCGE_GetIntFromCommandLine   (const char *para, GCGE_INT argc, char* argv[]);
GCGE_DOUBLE GCGE_GetDoubleFromCommandLine(const char *para, GCGE_INT argc, char* argv[]);
char *GCGE_GetCharFromCommandLine  (const char *para, GCGE_INT argc, char* argv[]);
void GCGE_Printf(char *fmt, ...);
/**
 * @brief 
 *
 * @param argc
 * @param argv
 * @param n
 * @param name
 */

#endif
