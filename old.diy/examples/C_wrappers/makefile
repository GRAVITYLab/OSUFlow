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
include ../../config_defs

INCLUDE = -I../../include
LIBS = -L../../lib -ldiy
OBJS = main.o
CCFLAGS =

.SUFFIXES: .c

.c.o:
	$(CC) -c $(CCFLAGS) $(INCLUDE) $<

default: all

all: example

example: $(OBJS)
	$(CC) -o example $(OBJS) $(LIBS) -lstdc++

clean:
	rm -f *.o example
