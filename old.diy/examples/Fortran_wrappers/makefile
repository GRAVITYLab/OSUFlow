#----------------------------------------------------------------------------
#
# Tom Peterka
# Argonne National Laboratory
# 9700 S. Cass Ave.
# Argonne, IL 60439
# tpeterka@mcs.anl.gov
#
# (C) 2011 by Argonne National Laboratory.
# See COPYRIGHT in top-level directory.
#
#----------------------------------------------------------------------------
INCLUDE = -I../../include
LIBS = -L../../lib -ldiy
OBJS = main.o
CC = mpif90
CCFLAGS =

.SUFFIXES: .f90

.f90.o:
	$(CC) -c $(CCFLAGS) $(INCLUDE) $<

default: all

all: example

example: $(OBJS)
	$(CC) -o example $(OBJS) $(LIBS) -lstdc++

clean:
	rm -f example *.o 
