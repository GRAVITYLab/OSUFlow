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

INCLUDE += -I. -I../include -I../src/bil-0.6.0/src
CCFLAGS = -c
CCFLAGS += -g

# architecture
ifeq ($(ARCH), LINUX)
CCFLAGS += -DLINUX
endif
ifeq ($(ARCH), BGP)
CCFLAGS += -DBGP
endif
ifeq ($(ARCH), BGQ)
CCFLAGS += -DBGQ
endif
ifeq ($(ARCH), XKE)
CCFLAGS += -DXKE
endif
ifeq ($(ARCH), MAC)
CCFLAGS += -DMAC
endif

# pnetcdf
ifeq ($(PNETCDF), NO)
CCFLAGS += -DDISABLE_PNETCDF
endif

# position-independent code
ifeq ($(FPIC), YES)
CCFLAGS += -fPIC
endif

# bil
OBJS = src/bil.o src/bil_misc.o src/bil_sched.o \
	src/bil_block.o src/bil_pio.o src/bil_timing.o

.SUFFIXES: .c

default: all

all: libbil.a

.c.o:
	$(CC) $(CCFLAGS) $(INCLUDE) $<

libbil.a : $(OBJS)
	rm -f $@
	ar cru $@ *.o

clean:
	rm -f *.o *.a
