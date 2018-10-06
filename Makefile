GCGEHOME = .
include $(GCGEHOME)/config/make.inc
.PHONY: clean


###########################################################
all: 	help


###########################################################
libs: 
	@echo "======================================="
	@echo "        Making libraries               "
	@cd $(GCGESRC);             $(MAKE) lib
	@cd $(GCGESRC)/app/csr;     $(MAKE) lib
	@cd $(EXTERNAL)/csr/src;    $(MAKE) lib
	@cd $(GCGESRC)/app/hypre;   $(MAKE) lib
#	@cd $(GCGESRC)/app/pase;    $(MAKE) lib
#	@cd $(GCGESRC)/app/petsc;   $(MAKE) lib
#	@cd $(GCGESRC)/app/slepc;   $(MAKE) lib

libgcge:
	@echo "======================================="
	@echo "        Making library GCGE             "
	@cd $(GCGEHOME)/src;         $(MAKE) lib

test-csr: 
	@cd $(GCGEHOME)/test/csr;    $(MAKE) exe run-solver

#test-hypre: 
#	@cd $(GCGEHOME)/test/hypre;  $(MAKE) exe run-ex1

#test-pase: 
#	@cd $(GCGEHOME)/test/pase;   $(MAKE) exe run-ex1

#test-petsc: 
#	@cd $(GCGEHOME)/test/petsc;  $(MAKE) exe run-ex1

#test-slepc: 
#	@cd $(GCGEHOME)/test/slepc;  $(MAKE) exe run-ex1

clean:
	@cd $(GCGEHOME)/src; make clean
	@cd $(GCGEHOME)/src/app/csr; make clean
	@cd $(GCGEHOME)/external/csr/src; make clean
	@cd $(GCGEHOME)/test/csr; make clean


cleanlibs: 
	@cd $(GCGELIB);  $(RM) $(RMFLAGS) *.a

help:
	@echo " "
	@echo "make {libs|clean|cleanlibs|help}" 
	@echo " "
	@echo "   libs         - create all libraries"
	@echo "   libgcge      - create GCGE library"
	@echo " "
	@echo "   test-csr     - test library GCGE      (only 1 test)"
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
