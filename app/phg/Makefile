########################################################################
GCGEHOME = ../..
include $(GCGEHOME)/config/make.inc
.PHONY: all lib clean help

########################################################################
FILES   = gcge_app_phg.c
HEADERS = gcge_app_phg.h
OBJS    = ${FILES:.c=.o}

INC     = -I$(GCGEINC) -I$(PHGINC) $(PHGFLAGS)
LIB     = $(LIBGCGEPHG)

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
	@echo 'Creating header files $(GCGEINC)/$(HEADERS)'
	@$(ARCH) $(ARCHFLAGS) $(LIB) $(OBJS)
	@$(RANLIB) $(LIB)
	@$(CP) -fR $(HEADERS) $(GCGEINC)

$(OBJS): $(HEADERS)
########################################################################
