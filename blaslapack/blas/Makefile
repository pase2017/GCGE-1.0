########################################################################
GCGEHOME = ../..
include $(GCGEHOME)/config/make.inc
.PHONY: all lib clean help

########################################################################
FILES   = daxpy.f dcopy.f ddot.f  dgemm.f  dgemv.f dger.f  dnrm2.f drot.f  drotg.f dscal.f dswap.f \
	  dsymm.f dsymv.f dsyr2.f dsyr2k.f dsyr.f  dsyrk.f dtrmm.f dtrmv.f dtrsm.f dtrsv.f idamax.f \
	  lsame.f xerbla.f
OBJS    = ${FILES:.f=.o}

INC     = 
LIB     = $(LIBBLAS)

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
	@$(ARCH) $(ARCHFLAGS) $(LIB) $(OBJS)
	@$(RANLIB) $(LIB)

########################################################################
