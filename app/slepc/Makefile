########################################################################
GCGEHOME = ../..
include $(GCGEHOME)/config/make.inc

########################################################################
OBJS = gcge_app_slepc.o 
LIBS = $(LIBSLEPC)
INC  = -I$(APPINC) -I$(GCGEINC) -I$(SLEPCINC)
LIB  = $(LIBGCGESLEPC)

########################################################################
all:	help

lib:	$(LIB) info

clean:
	@$(RM) $(RMFLAGS) *.o

help:
	@echo "make {lib|clean|help}"
	@echo " "
	@echo "   lib   - create library libgcgeslepc-1.X.a"
	@echo " "
	@echo "   clean - remove temporary files"
	@echo "   help  - print this message"
	@echo " "

info:
	@echo " "
	@echo "Summary:"
	@echo "Executables are located in bin/"
	@echo "Libraries are located in lib/"
	@echo " "


########################################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INC)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	$(CC)  $(CFLAGS) -c $*.c -o $*.o $(INC)
.cpp.o:
	@echo '$(CPP)  -c $(CPPFLAGS) ' $*.cpp
	@$(CPP)  $(CPPFLAGS) -c $*.cpp -o $*.o $(INC)

########################################################################
$(LIB): $(OBJS)
	@echo ''
	@echo 'Creating library $(LIB)'
	@$(ARCH) $(ARCHFLAGS) $(LIB)  $(OBJS)
	@$(RANLIB) $(LIB)

########################################################################
