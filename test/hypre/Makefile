############################################################
GCGEHOME = ../..
include $(GCGEHOME)/config/make.inc

#CC      = $(MPICC)
#CLINKER = $(MPICC)
############################################################
OBJS1 = test_matvec.o
OBJS2 = test_hypre_solver.o
OBJS3 = test_hypre_solver_ls.o
OBJS4 = test_gcge_lobpcg.o
OBJS5 = test_bcg.o

EXE1  = $(GCGEBIN)/test_matvec.exe 
EXE2  = $(GCGEBIN)/test_hypre_solver.exe
#EXE3  = $(GCGEBIN)/test_hypre_solver_ls.exe
EXE4  = $(GCGEBIN)/test_gcge_lobpcg.exe
EXE5  = $(GCGEBIN)/test_bcg.exe

# use DHAVE_CONFIG_H and DHYPRE_TIMING to use hypre_time
INC   = -I$(APPINC) -I$(GCGEINC) -I$(HYPREINC) -DHAVE_CONFIG_H -DHYPRE_TIMING
LIBS  = $(LIBGCGEHYPRE) $(LIBGCGE) $(LIBHYPRE) $(LIBLAPACK) $(LIBBLAS) $(LIBGFORTRAN) -lm

############################################################

############################################################

all:	    help

exe:         $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) 

clean:
	@$(RM) $(RMFLAGS) *.o $(GCGEBIN)/*.exe *.~

help:
	@echo " "
	@echo "make {exe|run-mesh-ani2D|clean|help}"
	@echo " "
	@echo "   exe   - compile two executables"
	@echo " "
	@echo "   xx    - test xx"
	@echo " "
	@echo "   clean - remove temporary files"
	@echo "   help  - print this message"
	@echo " "

info:
	@echo " "
	@echo "Summary:"
	@echo "Executables are located in bin/"
	@echo "Libraries   are located in lib/"
	@echo " "

run-matvec:
	$(DEBUGEXEC) $(DEBUGFLAGS) $(EXE1)

run-hypre-solver:
	$(DEBUGEXEC) $(DEBUGFLAGS) $(EXE2)

run-hypre-solver-ls:
	$(DEBUGEXEC) $(DEBUGFLAGS) $(EXE3)

run-gcge-lobpcg:
	$(DEBUGEXEC) $(DEBUGFLAGS) $(EXE4)

run-bcg:
	$(DEBUGEXEC) $(DEBUGFLAGS) $(EXE5)

############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INC)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	$(CC)  $(CFLAGS) -c $*.c -o $*.o $(INC)
.cpp.o:
	@echo '$(CPP)  -c $(CPPFLAGS) ' $*.cpp
	@$(CPP)  $(CPPFLAGS) -c $*.cpp -o $*.o $(INC)

############################################################
$(EXE1): $(OBJS1)
	@echo ''
	@echo 'Creating executable bin/xx.exe'
	$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS)
	@echo ''

$(EXE2): $(OBJS2)
	@echo ''
	@echo 'Creating executable bin/xxx.exe'
	$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS)
	@echo ''

$(EXE3): $(OBJS3)
	@echo ''
	@echo 'Creating executable bin/xxx.exe'
	$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS)
	@echo ''

$(EXE4): $(OBJS4) laplace.o 
	@echo ''
	@echo 'Creating executable bin/xxx.exe'
	$(CLINKER) $(LDFLAGS) -o $@ $^ $(LIBS)
	@echo ''

$(EXE5): $(OBJS5) laplace.o 
	@echo ''
	@echo 'Creating executable bin/xxx.exe'
	$(CLINKER) $(LDFLAGS) -o $@ $^ $(LIBS)
	@echo ''

############################################################
