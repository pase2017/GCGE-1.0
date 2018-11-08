GCGEHOME = .
include $(GCGEHOME)/config/make.inc
.PHONY: all libs clean help

###########################################################
all: 	help

###########################################################
libs: libgcge_core libgcge_csr libgcge_hypre

libgcge_core:
	@echo "======================================="
	@echo "        Making library GCGE CORE       "
	@cd $(GCGEHOME)/src;         $(MAKE) lib

libgcge_csr:
	@echo "======================================="
	@echo "        Making library GCGE             "
	@cd $(GCGEHOME)/csr;        $(MAKE) lib

libgcge_hypre:
	@echo "======================================="
	@echo "        Making library HYPRE           "
	@cd $(GCGEHOME)/app/hypre;  $(MAKE) lib

libgcge_petsc:
	@echo "======================================="
	@echo "        Making library PETSC           "
	@cd $(GCGEHOME)/app/petsc;  $(MAKE) lib

libgcge_slepc:
	@echo "======================================="
	@echo "        Making library SLEPC           "
	@cd $(GCGEHOME)/app/slepc;  $(MAKE) lib


#test-csr: 
#	@cd $(GCGEHOME)/test/csr;    $(MAKE) exe run-solver

#test-hypre: 
#	@cd $(GCGEHOME)/test/hypre;  $(MAKE) exe run-ex1

#test-pase: 
#	@cd $(GCGEHOME)/test/pase;   $(MAKE) exe run-ex1

#test-petsc: 
#	@cd $(GCGEHOME)/test/petsc;  $(MAKE) exe run-ex1

#test-slepc: 
#	@cd $(GCGEHOME)/test/slepc;  $(MAKE) exe run-ex1

clean:
	@cd $(GCGESRC);             make clean; rm -rf *~
	@cd $(GCGEHOME)/csr;        make clean; rm -rf *~
	@cd $(GCGEHOME)/app/hypre;  make clean; rm -rf *~
#	@cd $(GCGEHOME)/app/petsc;  make clean; rm -rf *~
#	@cd $(GCGEHOME)/app/slepc;  make clean; rm -rf *~
	@cd $(GCGEHOME)/test/hypre; make clean; rm -rf *~
#	@cd $(GCGEHOME)/test/slepc; make clean; rm -rf *~
#	@cd $(GCGEHOME)/test/petsc; make clean; rm -rf *~


cleanall: clean
	@cd $(GCGELIB);  $(RM) $(RMFLAGS) *.a
	@cd $(GCGEINC);  $(RM) $(RMFLAGS) *.h

help:
	@echo " "
	@echo "make {libs|clean|cleanlibs|help}" 
	@echo " "
	@echo "   libs         - create all libraries"
	@echo "   libgcge      - create GCGE library"
	@echo " "
#	@echo "   test-csr     - test library GCGE      (only 1 test)"
#	@echo "   test-hypre   - test library GCGEHypre (only 1 test)"
#	@echo "   test-pase    - test library GCGEPASE  (only 1 test)"
#	@echo "   test-petsc   - test library GCGEPetsC (only 1 test)"
#	@echo "   test-slepc   - test library GCGESlepC (only 1 test)"
	@echo " "
	@echo "   clean        - remove temporary files except libraries"
	@echo "   cleanlibs    - remove libraries"
	@echo " "
	@echo "   help         - print this information"
	@echo " "
