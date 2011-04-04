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

INCLUDE += -I../../include -I../../src/bil-0.6.0/src
LIBS += -L../../lib -ldiy -L../../src/bil-0.6.0 -lbil
OBJS = main.o
CCFLAGS =

.SUFFIXES: .cpp

.cpp.o:
	$(CXX) -c $(CCFLAGS) $(INCLUDE) $<

default: all

all: example

example: $(OBJS)
	$(CXX) -o example $(OBJS) $(LIBS)

clean:
	rm -f *.o 
