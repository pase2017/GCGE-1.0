############################################################
GCGEHOME = ../..
include $(GCGEHOME)/config/make.inc

############################################################
OBJS1 = test_solver.o
OBJS2 = test_matvec.o
OBJS3 = test_orthogonal.o
OBJS4 = test_xpw.o
OBJS5 = test_ops.o
OBJS6 = test_cg.o
OBJS7 = test_rayleighritz.o
OBJS8 = test_solver_dsyev.o
OBJS9 = test_orth_ops.o
OBJS10 = test_orth_subspace.o
OBJS11 = test_p.o
OBJS12 = test_w.o
EXE1  = $(GCGEBIN)/solver.exe 
EXE2  = $(GCGEBIN)/test_matvec.exe 
EXE3  = $(GCGEBIN)/test_orthogonal.exe 
EXE4  = $(GCGEBIN)/test_xpw.exe 
EXE5  = $(GCGEBIN)/test_ops.exe 
EXE6  = $(GCGEBIN)/test_cg.exe 
EXE7  = $(GCGEBIN)/test_rayleighritz.exe
EXE8  = $(GCGEBIN)/test_solver_dsyev.exe
EXE9  = $(GCGEBIN)/test_orth_ops.exe
EXE10  = $(GCGEBIN)/test_orth_subspace.exe
EXE11  = $(GCGEBIN)/test_p.exe
EXE12  = $(GCGEBIN)/test_w.exe

INC   = -I$(APPINC) -I$(GCGEINC) -I$(CSRINC)
LIBS  = $(LIBGCGECSR) $(LIBGCGE) $(LIBCSR) $(LIBLAPACK) \
	$(LIBBLAS) $(LIBGFORTRAN) -lm 

############################################################
MAT_A = $(GCGEHOME)/test/data/testA
MAT_M = $(GCGEHOME)/test/data/testB
#PARA1 = -file_A $(MAT_A) -file_B $(MAT_M) -gcge_ev_max_it 2 -gcge_nev 1
PARA = -file_A $(MAT_A) -file_B $(MAT_M)

############################################################

all:	    help

#exe:        $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9)
exe:        $(EXE1)

clean:
	@$(RM) $(RMFLAGS) *.o $(GCGEBIN)/*.exe 

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
	
run-matvec:
	$(DEBUGEXEC) $(DEBUGFLAGS) $(EXE2)

run-orth:
	$(DEBUGEXEC) $(DEBUGFLAGS) $(EXE3)

run-xpw:
	$(DEBUGEXEC) $(DEBUGFLAGS) $(EXE4)

run-ops:
	$(DEBUGEXEC) $(DEBUGFLAGS) $(EXE5)

run-cg:
	$(DEBUGEXEC) $(DEBUGFLAGS) $(EXE6)

run-rr:
	$(DEBUGEXEC) $(DEBUGFLAGS) $(EXE7)

run-dsyev:
	$(DEBUGEXEC) $(DEBUGFLAGS) $(EXE8)

run-orth-ops:
	$(EXE9)

run-orth-sub:
	$(EXE10)

run-p:
	valgrind -v $(EXE11)

run-w:
	valgrind -v $(EXE12)
	 
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

$(EXE4): $(OBJS4)
	@echo ''
	@echo 'Creating executable bin/xxx.exe'
	$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS)
	@echo ''

$(EXE5): $(OBJS5)
	@echo ''
	@echo 'Creating executable bin/xxx.exe'
	$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS)
	@echo ''

$(EXE6): $(OBJS6)
	@echo ''
	@echo 'Creating executable bin/xxx.exe'
	$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS)
	@echo ''

$(EXE7): $(OBJS7)
	@echo ''
	@echo 'Creating executable bin/xxx.exe'
	$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS)
	@echo ''

$(EXE8): $(OBJS8)
	@echo ''
	@echo 'Creating executable bin/xxx.exe'
	$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS)
	@echo ''

$(EXE9): $(OBJS9)
	@echo ''
	@echo 'Creating executable bin/xxx.exe'
	$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS)
	@echo ''

$(EXE10): $(OBJS10)
	@echo ''
	@echo 'Creating executable bin/xxx.exe'
	$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS)
	@echo ''

$(EXE11): $(OBJS11)
	@echo ''
	@echo 'Creating executable bin/xxx.exe'
	$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS)
	@echo ''

$(EXE12): $(OBJS12)
	@echo ''
	@echo 'Creating executable bin/xxx.exe'
	$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS)
	@echo ''


############################################################
