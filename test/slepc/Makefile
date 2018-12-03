############################################################
GCGEHOME = ../..
include $(GCGEHOME)/config/make.inc

MAIN2 = test_orthogonal.c
OBJS2 = test_orthogonal.o
EXE2 = bin_test_orthogonal.exe

MAIN1 = test_solver.c
OBJS1 = test_solver.o
EXE1 = ../../bin/test_solver.exe

MAIN1 = test_petsc_cg.c
OBJS1 = test_petsc_cg.o
EXE1 = ../../bin/test_petsc_cg.exe

MAIN1 = test_pcg.c
OBJS1 = test_pcg.o
EXE1 = ../../bin/test_pcg.exe

MAIN1 = test_bcg.c
OBJS1 = test_bcg.o
EXE1 = test_bcg.exe

MAIN1 = test_slepc_solver.c
OBJS1 = test_slepc_solver.o
EXE1 = ./test_slepc_solver.exe

INC   = -I$(GCGEINC) -I$(SLEPCINC) -I$(PETSCINC)
LIBS  = $(LIBGCGESLEPC) $(LIBGCGE) $(LIBSLEPC)

all: $(OBJS1) $(EXE1)
$(EXE1):$(OBJS1)
	$(CLINKER) -o $(EXE1) $(OBJS1) $(LIBS)

$(OBJS1):$(MAIN1)
	$(CC) $(CFLAGS) -c $(MAIN1) -o $(OBJS1) $(INC)

run: $(EXE1)
	$(EXEC) $(EXECFLAGS) $(EXE1) -mat_A ../data/A_5.petsc.bin -mat_B ../data/M_5.petsc.bin
#include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
