############################################################
GCGEHOME = ../..
include ${PETSC_DIR}/lib/petsc/conf/variables
include $(GCGEHOME)/config/make.inc

############################################################
OBJS1 = test_petsc_solver.o

EXE1  = $(GCGEBIN)/test_petsc_solver.exe 

INC   = -I$(APPINC) -I$(GCGEINC) -I$(PETSCINC)
LIBS  = $(LIBGCGEPETSC) $(LIBGCGE) $(LIBPETSC) $(LIBLAPACK) \
	$(LIBBLAS) $(LIBGFORTRAN) -lm 

############################################################
MAT_A = $(GCGEHOME)/test/data/testA
MAT_M = $(GCGEHOME)/test/data/testB
#PARA1 = -file_A $(MAT_A) -file_B $(MAT_M) -gcge_ev_max_it 2 -gcge_nev 1
PARA = -file_A $(MAT_A) -file_B $(MAT_M)

############################################################

all:	    help

exe:         $(EXE1)

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

run-solver:
	$(EXE1)
	
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

############################################################
