GCGEHOME = .
include $(GCGEHOME)/config/make.inc
.PHONY: all libs clean help

WITH_MPI = $(findstring mpi, $(CC))

install_libs =
install_libs += $(if $(WITHMPI) without-mpi  ,with-mpi  ) 
install_libs += libgcge_core libgcge_csr
install_libs += $(if $(HYPREINC) ,libgcge_hypre) 
install_libs += $(if $(PETSCINC) ,libgcge_petsc) 
install_libs += $(if $(SLEPCINC) ,libgcge_slepc) 

###########################################################
all: 	help

###########################################################

libs: $(install_libs)
	@echo $(install_libs)
	@echo "======================================="
	@echo "Build $(strip $(subst -, , $(install_libs))) complete."
	@echo "======================================="
	@echo "Now to install the library to $(INSTALLDIR) do:"
	@echo "make install"
	
without-mpi:
	@sed -i "s/#define GCGE_USE_MPI [01]/#define GCGE_USE_MPI 0/" $(GCGESRC)/gcge_type.h

with-mpi:
	@sed -i "s/#define GCGE_USE_MPI [01]/#define GCGE_USE_MPI 1/" $(GCGESRC)/gcge_type.h

libgcge_core:
	@echo "======================================="
	@echo "        Making library GCGE CORE       "
	@cd $(GCGEHOME)/src;        $(MAKE) lib

libgcge_csr:
	@echo "======================================="
	@echo "        Making library CSR             "
	@cd $(GCGEHOME)/app/csr;    $(MAKE) lib

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
	@cd $(GCGEHOME)/app/csr;    make clean; rm -rf *~
	@cd $(GCGEHOME)/app/hypre;  make clean; rm -rf *~
	@cd $(GCGEHOME)/app/petsc;  make clean; rm -rf *~
	@cd $(GCGEHOME)/app/slepc;  make clean; rm -rf *~
#	@cd $(GCGEHOME)/test/hypre; make clean; rm -rf *~
#	@cd $(GCGEHOME)/test/slepc; make clean; rm -rf *~
#	@cd $(GCGEHOME)/test/petsc; make clean; rm -rf *~


cleanall: clean
	@cd $(GCGELIB);  $(RM) $(RMFLAGS) *.a
	@cd $(GCGEINC);  $(RM) $(RMFLAGS) *.h

install:
	@echo "Create $(INSTALLDIR)/include and $(INSTALLDIR)/lib"
	@mkdir -p  $(INSTALLDIR)
	@mkdir -p  $(INSTALLDIR)/include
	@mkdir -p  $(INSTALLDIR)/lib
	@$(CP) -fR $(GCGEINC)/*.h $(INSTALLDIR)/include
	@$(AR) -x  $(GCGELIB)/*.a
	@$(ARCH) $(ARCHFLAGS) libGCGE-$(version).a *.o
	@$(RANLIB) libGCGE-$(version).a
	@$(RM) $(RMFLAGS) *.o
	@$(CP) -fR libGCGE-$(version).a $(INSTALLDIR)/lib
	@echo "Install complete."
	@echo "======================================="

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
	@echo "   cleanall     - remove libraries"
	@echo " "
	@echo "   help         - print this information"
	@echo " "
