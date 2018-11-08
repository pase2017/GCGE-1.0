############################################################
INSTALLDIR   = ${HOME}/soft/GCGE-2.0

############################################################
# DEFINE libs and dirs of include of OEM package
############################################################

LIBLAPACK    = ${INTEL_DIR}/mkl/lib/intel64/libmkl_lapack95_lp64.a
LIBBLAS      = ${INTEL_DIR}/mkl/lib/intel64/libmkl_blas95_lp64.a

############################################################
# DEFINE libs and dirs of include of OEM package
############################################################

PETSCINC     = ${PETSC_DIR}/include
PETSCINC     = 
LIBPETSC     = ${PETSC_KSP_LIB}

SLEPCFLAGS   = 
SLEPCINC     = ${SLEPC_DIR}/include
LIBSLEPC     = ${SLEPC_SYS_LIB}

HYRPREFLAGS  = -DHAVE_CONFIG_H -DHYPRE_TIMING
HYPREINC     = ${HYPRE_DIR}/include
LIBHYPRE     = ${HYPRE_DIR}/lib/libHYPRE.a

############################################################
VALGRIND      =
VALGRINDFLAGS =
############################################################
# DEFINE compilers, linkers and execution
############################################################
F77          = mpiifort
CC           = mpiicc
CPP          = mpiicpc
EXEC         = mpiexec
EXECFLAGS    = -np 2

FLINKER      = $(F77)
CLINKER      = $(CC)
CPPLINKER    = $(CPP)

DEBUG        = -g
FFLAGS       = -O3 -fPIC $(DEBUG)
CFLAGS       = -O3 -fPIC $(DEBUG)
CPPFLAGS     = -O3 -fPIC $(DEBUG)

DEBUGEXEC    = $(VALGRIND)
DEBUGFLAGS   = $(VALGRINDFLAGS)

############################################################
# DEFINE for mixed code
############################################################
# linking Fortran
LIBF         = -lgfortranbegin -lgfortran -lm -shared-libgcc
# linking CPP
LIBCPP       = -lstdc++
# other exist libs
LIBOE        = 
# linker flags
LDFLAGS      = -mkl=cluster

############################################################
# DEFINE dirs
############################################################

GCGEBIN       = $(GCGEHOME)/bin
GCGEINC       = $(GCGEHOME)/include
GCGELIB       = $(GCGEHOME)/lib
GCGESRC       = $(GCGEHOME)/src

############################################################
version       = 2.0

############################################################
# DEFINE libs and dirs of include
############################################################

LIBGCGE       = $(GCGELIB)/libgcge_core-$(version).a

LIBGCGECSR    = $(GCGELIB)/libgcge_csr-$(version).a
LIBGCGEHYPRE  = $(GCGELIB)/libgcge_hypre-$(version).a
LIBGCGEPASE   = $(GCGELIB)/libgcge_pase-$(version).a
LIBGCGEPHG    = $(GCGELIB)/libgcge_phg-$(version).a
LIBGCGEPETSC  = $(GCGELIB)/libgcge_petsc-$(version).a
LIBGCGESLEPC  = $(GCGELIB)/libgcge_slepc-$(version).a

############################################################

MAKE         = make
CP           = /bin/cp
RM           = /bin/rm
RMFLAGS      = -f
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib
############################################################
