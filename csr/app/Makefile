########################################################################
GCGEHOME = ../..
include $(GCGEHOME)/config/make.inc
########################################################################
CSRINC      = ../include
LIBGCGECSR  = ../lib/libgcge_csr.a

HEADER = gcge_app_csr.h
OBJS   = gcge_app_csr.o
INC    = -I$(GCGEINC) -I$(CSRINC)
LIB    = $(LIBGCGECSR)

########################################################################
all:	help

lib:	$(LIB)

clean:
	@$(RM) $(RMFLAGS) *.o

help:
	@echo "make {lib|clean|help}"
	@echo " "
	@echo "   lib   - create library $(LIB)"
	@echo " "
	@echo "   clean - remove temporary files"
	@echo "   help  - print this message"
	@echo " "

########################################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INC)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o $(INC)
.cpp.o:
	@echo '$(CPP)  -c $(CPPFLAGS) ' $*.cpp
	@$(CPP)  $(CPPFLAGS) -c $*.cpp -o $*.o $(INC)

########################################################################
$(LIB): $(OBJS)
	@echo ''
	@echo 'Creating library $(LIB)'
	@echo 'Creating library $(CSRINC)/$(HEADER)'
	@$(ARCH) $(ARCHFLAGS) $(LIB)  $(OBJS)
	@$(RANLIB) $(LIB)
	@$(CP) $(HEADER)  $(CSRINC)

########################################################################
