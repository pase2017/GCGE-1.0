/*
 * =====================================================================================
 *
 *       Filename:  gcge_utilities.c
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
#include <sys/resource.h>

#include "gcge_utilities.h"

//获取当前时刻时间
GCGE_DOUBLE GCGE_GetTime()
{
    struct rusage usage;
    GCGE_DOUBLE ret;
      
    if(getrusage(RUSAGE_SELF, &usage) == -1) 
        printf("Error in GCGE_GetTime!\n");
 
    ret = ((GCGE_DOUBLE) usage.ru_utime.tv_usec)/1000000;
 
    ret += usage.ru_utime.tv_sec;
 
    return ret;
}

//从命令行读取int型参数
GCGE_INT GCGE_GetIntFromCommandLine(const char *para, GCGE_INT argc, char* argv[])
{
    GCGE_INT arg_index = 0, intpara;
    while(arg_index < argc) 
    {
        if(0 == strcmp(argv[arg_index], para)) 
        {
            intpara = atoi(argv[++arg_index]);
            break;
        }
        arg_index++;
    }
    return intpara;
}

//从命令行读取double型参数
GCGE_DOUBLE GCGE_GetDoubleFromCommandLine(const char *para, GCGE_INT argc, char* argv[])
{
    GCGE_INT arg_index = 0;
    GCGE_DOUBLE doublepara;
    while(arg_index < argc) 
    {
        if(0 == strcmp(argv[arg_index], para)) 
        {
            doublepara = atof(argv[++arg_index]);
            break;
        }
        arg_index++;
    }
    return doublepara;
}

//从命令行读取char*型参数
char* GCGE_GetCharFromCommandLine(const char *para, GCGE_INT argc, char* argv[])
{
    GCGE_INT arg_index = 0;
    char *charpara;
    while(arg_index < argc) 
    {
        if(0 == strcmp(argv[arg_index], para)) 
        {
            //charpara = (const char *)(&(argv[++arg_index]));
            printf("argv: %s\n", argv[arg_index+1]);
            charpara = argv[++arg_index];
            printf("charpara: %s\n", charpara);
            break;
        }
        arg_index++;
    }
    return charpara;
}


void GCGE_Printf(char *fmt, ...)
{
#if GCGE_USE_MPI
    GCGE_INT myrank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if(0 == myrank) {
        va_list vp;
        va_start(vp, fmt);
        vprintf(fmt, vp);
        va_end(vp);
    }
#else
    va_list vp;
    va_start(vp, fmt);
    vprintf(fmt, vp);
    va_end(vp);
#endif
}

