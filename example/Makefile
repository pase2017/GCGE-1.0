############################################################
GCGEHOME = ..
include $(GCGEHOME)/config/make.inc
.PHONY: all clean help
############################################################

EXE_CSR    = csr_solver.exe
EXE_HYPRE  = hypre_solver.exe
EXE_PETSC  = petsc_solver.exe
EXE_PHG    = phg_solver.exe
EXE_SLEPC  = slepc_solver.exe
           
OBJ_CSR    = csr_solver.o
OBJ_HYPRE  = hypre_solver.o
OBJ_PETSC  = petsc_solver.o
OBJ_PHG    = phg_solver.o
OBJ_SLEPC  = slepc_solver.o
           
C_CSR      = csr_solver.c
C_HYPRE    = hypre_solver.c
C_PETSC    = petsc_solver.c
C_PHG      = phg_solver.c
C_SLEPC    = slepc_solver.c
           
INC_CSR    = $(CSRFLAGS)   -I$(GCGEINC)            
INC_HYPRE  = $(HYPREFLAGS) -I$(GCGEINC) -I$(HYPREINC) 
INC_PETSC  = $(PETSCFLAGS) -I$(GCGEINC) -I$(PETSCINC) 
INC_PHG    = $(PHGFLAGS)   -I$(GCGEINC) -I$(PHGINC)   
INC_SLEPC  = $(SLEPCFLAGS) -I$(GCGEINC) -I$(SLEPCINC) 
          
LIBS_CSR   = $(LIBGCGECSR)   $(LIBGCGE) $(LIBLAPACK) $(LIBBLAS) $(LIBF)
LIBS_HYPRE = $(LIBGCGEHYPRE) $(LIBGCGE) $(LIBHYPRE)
LIBS_PETSC = $(LIBGCGEPETSC) $(LIBGCGE) $(LIBPETSC)
LIBS_PHG   = $(LIBGCGEPHG)   $(LIBGCGE) $(LIBPHG)
LIBS_SLEPC = $(LIBGCGESLEPC) $(LIBGCGE) $(LIBSLEPC)

############################################################

all:	    help

clean:
	@$(RM) $(RMFLAGS) *.o *.exe

help:
	@echo " "
	@echo "make {run-xx-solver|clean|help}"
	@echo " "
	@echo "   run-csr-solver   - test gcge solver in csr"
	@echo "   run-hypre-solver - test gcge solver in hypre"
	@echo "   run-petsc-solver - test gcge solver in petsc"
	@echo "   run-phg-solver   - test gcge solver in phg"
	@echo "   run-slepc-solver - test gcge solver in slepc"
	@echo " "
	@echo "   clean            - remove temporary files"
	@echo "   help             - print this message"
	@echo " "

############################################################

run-csr-solver: $(EXE_CSR)
	$(DEBUGEXEC) $(DEBUGFLAGS) ./$(EXE_CSR)

run-hypre-solver: $(EXE_HYPRE)
	$(DEBUGEXEC) $(DEBUGFLAGS) ./$(EXE_HYPRE)

run-petsc-solver: $(EXE_PHG)
	$(DEBUGEXEC) $(DEBUGFLAGS) ./$(EXE_PETSC)

run-phg-solver: $(EXE_PHG)
	$(DEBUGEXEC) $(DEBUGFLAGS) ./$(EXE_PHG)

run-slepc-solver: $(EXE_SLEPC)
	$(DEBUGEXEC) $(DEBUGFLAGS) ./$(EXE_SLEPC)

############################################################

$(EXE_CSR): $(OBJ_CSR)
	@echo ''
	@echo 'Creating executable $(EXE_CSR)'
	@$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS_CSR)
	@echo ''

$(EXE_HYPRE): $(OBJ_HYPRE)
	@echo ''
	@echo 'Creating executable $(EXE_HYPRE)'
	@$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS_HYPRE)
	@echo ''

$(EXE_PETSC): $(OBJ_PETSC)
	@echo ''
	@echo 'Creating executable $(EXE_PETSC)'
	@$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS_PETSC)
	@echo ''

$(EXE_PHG): $(OBJ_PHG)
	@echo ''
	@echo 'Creating executable $(EXE_PHG)'
	@$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS_PHG)
	@echo ''

$(EXE_SLEPC): $(OBJ_SLEPC)
	@echo ''
	@echo 'Creating executable $(EXE_SLEPC)'
	@$(CLINKER) $(LDFLAGS) -o $@ $< $(LIBS_SLEPC)
	@echo ''

############################################################

$(OBJ_CSR):$(C_CSR)
	@echo '$(CC) $(CFLAGS) -c $(C_CSR)'
	@$(CC)  $(CFLAGS) -c $(C_CSR) -o $(OBJ_CSR) $(INC_CSR)

$(OBJ_HYPRE):$(C_HYPRE)
	@echo '$(CC) $(CFLAGS) -c $(C_HYPRE)'
	@$(CC)  $(CFLAGS) -c $(C_HYPRE) -o $(OBJ_HYPRE) $(INC_HYPRE)

$(OBJ_PETSC):$(C_PETSC)
	@echo '$(CC) $(CFLAGS) -c $(C_PETSC)'
	@$(CC)  $(CFLAGS) -c $(C_PETSC) -o $(OBJ_PETSC) $(INC_PETSC)

$(OBJ_PHG):$(C_PHG)
	@echo '$(CC) $(CFLAGS) -c $(C_PHG)'
	@$(CC)  $(CFLAGS) -c $(C_PHG) -o $(OBJ_PHG) $(INC_PHG)

$(OBJ_SLEPC):$(C_SLEPC)
	@echo '$(CC) $(CFLAGS) -c $(C_SLEPC)'
	@$(CC)  $(CFLAGS) -c $(C_SLEPC) -o $(OBJ_SLEPC) $(INC_SLEPC)

############################################################
