############################################################
GCGEHOME = ../..
include $(GCGEHOME)/config/make.inc

############################################################
OBJS2 = test_phg_solver.o

EXE2  = test_phg_solver.exe

INC   = -I$(GCGEINC) -I$(PHGINC) $(PHGFLAGS)
LIBS  = $(LIBGCGEPHG) $(LIBGCGE) $(LIBPHG)
############################################################

############################################################

all:	    help

#exe:        $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) 
exe:        $(EXE2)

clean:
	@$(RM) $(RMFLAGS) *.o *.exe

help:
	@echo " "
	@echo "make {exe|run-mesh-ani2D|clean|help}"
	@echo " "
	@echo "   exe   - compile all executables"
	@echo " "
	@echo "   xx    - test xx"
	@echo " "
	@echo "   clean - remove temporary files"
	@echo "   help  - print this message"
	@echo " "

run-phg-solver:
	$(DEBUGEXEC) $(DEBUGFLAGS) $(EXE2)

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
$(EXE2): $(OBJS2)
	@echo ''
	@echo 'Creating executable bin/xxx.exe'
	$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS)
	@echo ''

############################################################
