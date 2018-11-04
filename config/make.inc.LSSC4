############################################################
# DEFINE libs and dirs of include of OEM package
############################################################

MKLPATH      = /share/soft/intel_2018_update1/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64
LIBLAPACK    = $(MKLPATH)/libmkl_lapack95_lp64.a
LIBBLAS      = $(MKLPATH)/libmkl_blas95_lp64.a
LIBGFORTRAN  = 

############################################################
# DEFINE libs and dirs of include of OEM package
############################################################

CSRINC       = $(EXTERNAL)/csr/include
LIBCSR       = $(EXTERNAL)/csr/lib/libcsr.a -lm

PETSCINC     = $(EXTERNAL)/petsc/include
LIBPETSC     = ${PETSC_KSP_LIB} $(EXTERNAL)/petsc/lib/libpetsc.a

SLEPCINC     = 
LIBSLEPC     = 

HYPREINC     = /share/home/hhxie/liyu/hypre-2.11.2/src/hypre/include
LIBHYPRE     = /share/home/hhxie/liyu/hypre-2.11.2/src/hypre/lib/libHYPRE.a

EXTERNALINCLUDEDIRS = -I$(PETSCINC) -I$(SLEPCINC) -I$(HYPREINC)
EXTERNALLINKLIBS    = $(LIBPETSC) $(LIBSLEPC) $(LIBHYPRE)

############################################################
MPIF77       = mpiifort
MPICC        = mpiicc
MPICPP       = mpiicpc
MPIEXEC      = mpiexec
MPIEXECFLAGS = -np 2
############################################################
VALGRIND       = 
VALGRINDFLAGS  =
VALGRINDFLAGSN =
############################################################
# DEFINE compilers, linkers and execution
############################################################
F77          = $(MPIF77)
CC           = $(MPICC)
CPP          = $(MPICPP)
MAKE         = make

FLINKER      = $(F77) -mkl=cluster
CLINKER      = $(CC)  -mkl=cluster
CPPLINKER    = $(CPP) -mkl=cluster

DEBUG        = -g
FFLAGS       = -O2       $(DEBUG)
CFLAGS       = -O2 -fPIC $(DEBUG)
CPPFLAGS     = -O2 -fPIC $(DEBUG)

DEBUGEXEC    = $(VALGRIND)
DEBUGFLAGS   = $(VALGRINDFLAGS)

EXEC         = 
EXECFLAGS    = 

############################################################
# DEFINE for mixed code
############################################################
# linking Fortran
LIBF         = -lm
# linking CPP
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
