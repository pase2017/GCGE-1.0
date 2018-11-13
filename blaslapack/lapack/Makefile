########################################################################
GCGEHOME = ../..
include $(GCGEHOME)/config/make.inc
.PHONY: all lib clean help

########################################################################
FILES   = dgesv.f  dgetf2.f dgetrf.f   dgetrs.f disnan.f dlacpy.f dlae2.f  dlaebz.f dlaev2.f \
	  dlagtf.f dlagts.f dlaisnan.f dlanst.f dlansy.f dlapy2.f dlarfb.f dlarf.f  dlarfg.f \
	  dlarft.f dlarnv.f dlartg.f   dlaruv.f dlascl.f dlaset.f dlasr.f  dlasrt.f dlassq.f \
	  dlaswp.f dlasyf.f dlatrd.f   dorg2l.f dorg2r.f dorgql.f dorgqr.f dorgtr.f dorm2l.f \
	  dorm2r.f dormql.f dormqr.f   dormtr.f dposv.f  dpotf2.f dpotrf.f dpotrs.f dstebz.f \
	  dstein.f dsteqr.f dsterf.f   dsyev.f  dsyevx.f dsysv.f  dsytd2.f dsytf2.f dsytrd.f \
	  dsytrf.f dsytrs.f ieeeck.f   iladlc.f iladlr.f ilaenv.f iparmq.f \
	  dlamch.f
OBJS    = ${FILES:.f=.o}

INC     = 
LIB     = $(LIBLAPACK)

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
