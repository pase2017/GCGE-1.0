############################################################
# DEFINE libs and dirs of include of OEM package
############################################################

#LIBLAPACK    = /home/nzhang/GCGE_newest/lib/liblapack.a
#LIBBLAS      = /home/nzhang/GCGE_newest/lib/blas_LINUX.a
LIBLAPACK    = -llapack -lgfortran -lm
LIBBLAS      = -lblas

############################################################
# DEFINE libs and dirs of include of OEM package
############################################################

CSRINC       = $(EXTERNAL)/csr/include
LIBCSR       = $(EXTERNAL)/csr/lib/libcsr.a -lm

PETSCINC     = 
LIBPETSC     = 

SLEPCINC     = 
LIBSLEPC     = 

HYPREINC     = 
LIBHYPRE     = 

EXTERNALINCLUDEDIRS = -I$(PETSCINC) -I$(SLEPCINC) -I$(HYPREINC)
EXTERNALLINKLIBS    = $(LIBPETSC) $(LIBSLEPC) $(LIBHYPRE) $(DEBUGLIBS)

############################################################
MPIF77       = mpif77
MPICC        = mpicc
MPICPP       = mpic++
MPIEXEC      = mpiexec
#MPIEXEC     = mpirun
MPIEXECFLAGS = -np 2
############################################################
VALGRIND      = valgrind
VALGRINDFLAGS = --tool=memcheck --leak-check=full -v --track-origins=yes --show-reachable=yes --read-var-info=yes
VALGRINDFLAGSN = --leak-check=full -v --show-reachable=yes 
############################################################
# DEFINE compilers, linkers and execution
############################################################
F77          = gfortran
CC           = gcc
CPP          = g++
MAKE         = make

FLINKER      = $(F77)
CLINKER      = $(CC)
CPPLINKER    = $(CPP)

DEBUG        = -g #-DMEMWATCH -DMW_STDIO
FFLAGS       = -O2 $(DEBUG)
CFLAGS       = -O2 -fPIC $(DEBUG)
CPPFLAGS     = -O2 -fPIC $(DEBUG)

DEBUGEXEC    = $(VALGRIND)
DEBUGFLAGS   = $(VALGRINDFLAGS)
#DEBUGLIBS    = -lmemwatch

EXEC         = 
EXECFLAGS    = 

############################################################
# DEFINE for mixed code
############################################################
#for linking Fortran
LIBF         = -lgfortranbegin -lgfortran -lm -shared-libgcc
#for linking CPP
LIBCPP       = -lstdc++
# other exist libs
LIBOE        = 
# lib dir
LDFLAGS      =

############################################################
# DEFINE dirs
############################################################

GCGEBIN       = $(GCGEHOME)/bin
GCGEINC       = $(GCGEHOME)/include
GCGELIB       = $(GCGEHOME)/lib
GCGESRC       = $(GCGEHOME)/src
GCGESRC       = $(GCGEHOME)/src

APPINC        = $(GCGEINC)/app

EXTERNAL      = $(GCGEHOME)/external
############################################################
version       = 1.0

############################################################
# DEFINE libs and dirs of include
############################################################

LIBGCGE       = $(GCGELIB)/libgcge_core-$(version).a
LIBGCGECSR    = $(GCGELIB)/libgcge_csr-$(version).a
LIBGCGEHYPRE  = $(GCGELIB)/libgcge_hypre-$(version).a
LIBGCGEPASE   = $(GCGELIB)/libgcge_pase-$(version).a
LIBGCGEPETSC  = $(GCGELIB)/libgcge_petsc-$(version).a
LIBGCGESLEPC  = $(GCGELIB)/libgcge_slepc-$(version).a

INCLUDE      = $(GCGEINC) $(APPINC)
############################################################

RM           = /bin/rm
RMFLAGS      = -f
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib
############################################################
