########################################################################
GCGEHOME = ..
include $(GCGEHOME)/config/make.inc
.PHONY: all lib clean help

########################################################################
FILES   = gcge_utilities.c gcge_statistics.c gcge_para.c gcge_ops.c \
          gcge_workspace.c \
          gcge_matvec.c \
          gcge_orthogonal.c gcge_rayleighritz.c gcge_cg.c gcge_xpw.c \
          gcge_eigsol.c \
          gcge_solver.c
HEADERS = gcge_utilities.h gcge_statistics.h gcge_para.h gcge_ops.h \
          gcge_workspace.h \
          gcge_matvec.h \
          gcge_orthogonal.h gcge_rayleighritz.h gcge_cg.h gcge_xpw.h \
          gcge_eigsol.h \
          gcge_solver.h \
          gcge.h
OBJS    = ${FILES:.c=.o}

INC     = 
LIB     = $(LIBGCGE)

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
	@$(CP) -fR gcge_*.h $(GCGEINC)
	@$(CP) -fR gcge.h   $(GCGEINC)

$(OBJS): $(HEADERS)
########################################################################
