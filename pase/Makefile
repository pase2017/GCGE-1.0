########################################################################
GCGEHOME = ..
include $(GCGEHOME)/config/make.inc
.PHONY: all lib clean help

########################################################################
FILES   = pase_matvec.c pase_mg.c pase_ops.c pase_solver.c
HEADERS = pase_matvec.h pase_mg.h pase_ops.h pase_solver.h
OBJS    = ${FILES:.c=.o}

INC     = 
LIB     = $(LIBPASE)

########################################################################
all:    help

lib:    $(LIB)

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
	@echo 'Creating header files $(GCGEINC)/{$(HEADERS)}'
	@$(ARCH) $(ARCHFLAGS) $(LIB) $(OBJS)
	@$(RANLIB) $(LIB)
#	@$(CP) -fR pase_*.h $(GCGEINC)
#	@$(CP) -fR pase.h   $(GCGEINC)

$(OBJS): $(HEADERS)
########################################################################
