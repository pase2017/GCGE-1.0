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

#include "gcge_type.h"

GCGE_DOUBLE GCGE_GetTime();
GCGE_INT    GCGE_GetIntFromCommandLine   (const char *para, GCGE_INT argc, char* argv[]);
GCGE_DOUBLE GCGE_GetDoubleFromCommandLine(const char *para, GCGE_INT argc, char* argv[]);
char *GCGE_GetCharFromCommandLine  (const char *para, GCGE_INT argc, char* argv[]);
/**
 * @brief 
 *
 * @param argc
 * @param argv
 * @param n
 * @param name
 */

#endif
