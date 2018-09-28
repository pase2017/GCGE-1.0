########################################################################
GCGEHOME = ..
include $(GCGEHOME)/config/make.inc

########################################################################
OBJS = gcge_utilities.o gcge_statistics.o gcge_para.o gcge_ops.o \
       gcge_matvec.o gcge_workspace.o \
	   gcge_orthogonal.o gcge_rayleighritz.o gcge_cg.o gcge_xpw.o \
	   gcge_eigsol.o \
	   gcge_solver.o
INC  = -I$(GCGEINC)
LIB  = $(LIBGCGE)
LIBS = $(LIBBLAS) $(LIBLAPACK)

########################################################################
all:	    help

lib:	    $(LIB) info

clean:
	    @$(RM) $(RMFLAGS) *.o $(GCGEBIN)/*.exe

help:
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib - create library libmesh-1.X.a"
	    @echo " "
	    @echo "   clean - clean the package"
	    @echo "   help  - print this message"
	    @echo " "

info:
	    @echo " "
	    @echo "Summary:"
	    @echo "Executables are located in bin/"
	    @echo "Libraries   are located in lib/"
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
