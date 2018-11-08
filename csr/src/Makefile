########################################################################
OBJS    = csr.o
HEADER  = csr.h
INC     = ../include
LIB     = ../lib/libcsr.a
CC      = gcc
CFLAG   = -g
########################################################################
CP           = cp
RM           = rm
RMFLAGS      = -f
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib
########################################################################
all:

lib:    $(LIB)

clean: 
	@$(RM) $(RMFLAGS) $(OBJS) $(LIB)

########################################################################
.c.o:
	@echo '$(CC) -c $(CFLAGS) ' $*.c
	@$(CC) $(CFLAGS) -c $*.c -o $*.o

########################################################################
$(LIB): $(OBJS)
	@echo ''
	@echo 'Creating library $(LIB)'
	@echo 'Creating header file $(INC)/$(HEADER)'
	@$(ARCH)   $(ARCHFLAGS) $(LIB)  $(OBJS)
	@$(RANLIB) $(LIB)
	@$(CP)     $(HEADER) $(INC)

########################################################################
