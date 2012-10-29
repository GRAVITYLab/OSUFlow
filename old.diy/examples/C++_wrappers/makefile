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
CXX = mpicxx
CCFLAGS =

.SUFFIXES: .cpp

.cpp.o:
	$(CXX) -c $(CCFLAGS) $(INCLUDE) $<

default: all

all: example

example: $(OBJS)
	$(CXX) -o example $(OBJS) $(LIBS)

clean:
	rm -f *.o example
