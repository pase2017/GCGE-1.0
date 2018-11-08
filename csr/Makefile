GCGEHOME = ../
include $(GCGEHOME)/config/make.inc
.PHONY: all lib clean help


###########################################################
all: 	help


###########################################################
lib: 
	@echo "======================================="
	@echo "        Making libraries               "
	@cd $(GCGEHOME)/csr/src;    $(MAKE) lib
	@cd $(GCGEHOME)/csr/app;    $(MAKE) lib

clean:
	@cd $(GCGEHOME)/csr/src;    $(MAKE) clean; rm -rf *~
	@cd $(GCGEHOME)/csr/app;    $(MAKE) clean; rm -rf *~


help:
	@echo " "
	@echo "make {libs|clean|cleanlibs|help}" 
	@echo " "
	@echo "   lib          - create all libraries"
	@echo " "
	@echo "   clean        - remove temporary files except libraries"
	@echo " "
	@echo "   help         - print this information"
	@echo " "
