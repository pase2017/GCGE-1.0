############################################################
INSTALLDIR   = ${HOME}/GCGE-2.0

############################################################
# DEFINE libs and dirs of include of OEM package
############################################################

LIBLAPACK    = $(GCGEHOME)/blaslapack/liblapack.a
LIBBLAS      = $(GCGEHOME)/blaslapack/libblas.a

############################################################
# DEFINE libs and dirs of include of OEM package
############################################################

PETSCFLAGS   = 
PETSCINC     = ${PETSC_INCLUDE}
LIBPETSC     = ${PETSC_KSP_LIB}

SLEPCFLAGS   = 
SLEPCINC     = ${SLEPC_INCLUDE}
LIBSLEPC     = ${SLEPC_SYS_LIB}

HYRPREFLAGS  = -DHAVE_CONFIG_H -DHYPRE_TIMING
HYPREINC     = /home/liyu/WorkSpace/hypre-2.11.2/src/hypre/include
LIBHYPRE     = /home/liyu/WorkSpace/hypre-2.11.2/src/hypre/lib/libHYPRE.a

############################################################
MPIF77       = mpif77
MPICC        = mpicc
MPICPP       = mpic++
MPIEXEC      = mpiexec
MPIEXECFLAGS = -np 2
############################################################
VALGRIND      = valgrind
VALGRINDFLAGS = --tool=memcheck --leak-check=full -v --track-origins=yes --show-reachable=yes --read-var-info=yes
############################################################
# DEFINE compilers, linkers and execution
############################################################
F77          = $(MPIF77)
CC           = $(MPICC)
CPP          = $(MPICPP)
MAKE         = make

FLINKER      = $(F77)
CLINKER      = $(CC)
CPPLINKER    = $(CPP)

DEBUG        = -g
FFLAGS       = -O2 -fPIC $(DEBUG)
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
LIBF         = -lgfortranbegin -lgfortran -lm -shared-libgcc
# linking CPP
LIBCPP       = -lstdc++
# other exist libs
LIBOE        = 
# linker flags
LDFLAGS      =

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

CP           = /bin/cp
RM           = /bin/rm
RMFLAGS      = -f
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib
############################################################
