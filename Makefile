#----------------------------------------------------------------------------
#
# makefile
#
# Copyright (c) 2009 Han-Wei Shen and Tom Peterka
#
# Contact:
#
# Han-Wei Shen
# The Ohio State University
# Columbus, OH
#
# Tom Peterka
# MCS Radix Lab
# Argonne National Laboratory
# 9700 S. Cass Ave.
# Argonne, IL 60439
# tpeterka@mcs.anl.gov
#
# All rights reserved. May not be used, modified, or copied
# without permission
#
#----------------------------------------------------------------------------

#ARCH = MAC_OSX
ARCH = LINUX
#ARCH = BGP
#ARCH = FD
#ARCH = EUREKA
#ARCH = BB

MPE = NO

LIBNAME = OSUFlow
RM = rm 
AR = ar cq

TOP = ..

### mac version ####

ifeq ($(ARCH),MAC_OSX)
C++ = g++
CC  = gcc
CCFLAGS = -g -c -DMAC_OSX -DGRAPHICS  -DDEBUG_MODE
LIBS  = -framework GLUT -framework OpenGL 
endif

### linux version ###

ifeq ($(ARCH),LINUX)
C++ = mpicxx
ifeq ($(MPE), YES)
C++   = mpecxx -mpilog
endif
THREADS = -fopenmp
CCFLAGS = -c -DLINUX -DMPI -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
CCFLAGS += -DGRAPHICS
CCFLAGS += -g
#CCFLAGS += -Wall -Wextra
LIBS = -lm -lglut -lGL
endif

### BG/P version ###

ifeq ($(ARCH), BGP)
INCLUDE = -I/usr/local/include \
	-I/usr/X11R6/include \
	-I/bgsys/drivers/ppcfloor/arch/include \

LIB  = -lm
C++ = mpixlcxx_r
ifeq ($(MPE), YES)
C++ = /home/chan/mpe_work/install_ibm/bin/mpecc -mpilog
endif
THREADS = -qsmp=omp:noauto
#THREADS = 
CCFLAGS += -O3 -qarch=450d -qtune=450
CCFLAGS += -c -DBGP -DMPI -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
ifeq ($(MPE), YES)
CCFLAGS += -DMPE
endif
endif

### sicortex version ###

ifeq ($(ARCH), FD)

INCLUDE = -I/usr/include
LIB  = -lm
C++   = mpicxx
ifeq ($(MPE), YES)
C++   = mpecxx -mpilog
endif
CCFLAGS = -c -DFD -DMPI -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
#CCFLAGS += -g3
CCFLAGS += -O3
#CCFLAGS += -Wall -Wextra
ifeq ($(MPE), YES)
CCFLAGS += -DMPE
endif
THREADS = -fopenmp

endif

### breadboard version ###

ifeq ($(ARCH), BB)

INCLUDE = -I/usr/include
LIB  = -lm
C++   = mpicxx
ifeq ($(MPE), YES)
C++   = mpecxx -mpilog
endif
CCFLAGS = -c -DBB -DMPI -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
#CCFLAGS += -g3
CCFLAGS += -O3
#CCFLAGS += -Wall -Wextra
ifeq ($(MPE), YES)
CCFLAGS += -DMPE
endif
THREADS = -fopenmp

endif

### eureka version ###

ifeq ($(ARCH), EUREKA)

INCLUDE = -I/usr/include
LIB  = -lm
C++   = mpicxx
ifeq ($(MPE), YES)
C++   = mpecxx -mpilog
endif
#CCFLAGS += -g3
CCFLAGS = -c -DEUREKA -DMPI -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX# -DGRAPHICS
CCFLAGS += -O3
#CCFLAGS += -Wall -Wextra
ifeq ($(MPE), YES)
CCFLAGS += -DMPE
endif
THREADS = -openmp

endif

##########

INCLUDE += -I.

OBJS =  Candidate.o  Grid.o  polynomials.o  TimeVaryingFieldLine.o \
	eigenvals.o  Interpolator.o  Rake.o	    Topology.o \
	eigenvecs.o  IsoSurf.o	     Solution.o     triangulator.o \
	Element.o    StreakLine.o    VectorMatrix.o \
	Field.o      PathLine.o      Streamline.o \
	FieldLine.o  Plot3DReader.o  TimeLine.o \
	OSUFlow.o    FileReader.o calc_subvolume.o \
	LatticeAMR.o Partition.o FlashAMR.o ComputeFieldLines.o \
	Lattice4D.o \

SRCS =  Candidate.C  Grid.C  polynomials.C  TimeVaryingFieldLine.C \
	eigenvals.C  Interpolator.C  Rake.C	    Topology.C \
	eigenvecs.C  IsoSurf.C	     Solution.C     triangulator.C \
	Element.C    StreakLine.C    VectorMatrix.C \
	Field.C      PathLine.C      Streamline.C \
	FieldLine.C  Plot3DReader.C  TimeLine.C \
	OSUFlow.C    FileReader.C calc_subvolume.C \
	LatticeAMR.C  Partition.C FlashAMR.C ComputeFieldLines.C \
	Lattice4D.C \

.SUFFIXES: .C

.C.o:
	$(C++) $(CCFLAGS) $(INCLUDE) $(THREADS)  $<

default: all

all: lib$(LIBNAME).a testmain testmain2 gldraw gldraw2 \
	testmainPathline testmainStreak gldrawPathline gldrawPathline2  gldrawStreak \
	gldrawStreak3 gldrawPathline3 gldrawPathline4 gldrawFlash gldrawFlash2 gldrawFlash3 \
	gldrawFlash4 gldrawFlashData gldrawFlashTime gldrawFlashPathline mpitest mpiamrtest

# deprecated - remove eventually
#	testmain3 testmain4 gldrawStreak2 gldraw3 gldraw4

lib$(LIBNAME).a : $(OBJS)
	$(RM) -f $@
	$(AR) $@ $(OBJS) 

testmain: testmain.o lib$(LIBNAME).a
	$(C++) -o testmain testmain.o -L. -l$(LIBNAME) -lm

testmain2: testmain2.o lib$(LIBNAME).a
	$(C++) -o testmain2 testmain2.o -L. -l$(LIBNAME) -lm

#testmain3: testmain3.o lib$(LIBNAME).a
#	$(C++) -o testmain3 testmain3.o -L. -l$(LIBNAME) -lm

#testmain4: testmain4.o lib$(LIBNAME).a
#	$(C++) -o testmain4 testmain4.o -L. -l$(LIBNAME) -lm

mpiamrtest: MpiAmrDraw.o lib$(LIBNAME).a
	$(C++) -o mpiamrtest MpiAmrDraw.o $(THREADS) -L. -l$(LIBNAME) $(LIBS) 

mpitest: MpiDraw.o lib$(LIBNAME).a
	$(C++) -o mpitest MpiDraw.o $(THREADS) -L. -l$(LIBNAME) $(LIBS) 

gldraw: gldraw.o  lib$(LIBNAME).a
	$(C++) -o gldraw gldraw.o -L. -l$(LIBNAME) $(LIBS) 

gldraw2: gldraw2.o  lib$(LIBNAME).a
	$(C++) -o gldraw2 gldraw2.o -L. -l$(LIBNAME) $(LIBS) 

#gldraw3: gldraw3.o  lib$(LIBNAME).a
#	$(C++) -o gldraw3 gldraw3.o -L. -l$(LIBNAME) $(LIBS)

#gldraw4: gldraw4.o  lib$(LIBNAME).a
#	$(C++) -o gldraw4 gldraw4.o -L. -l$(LIBNAME) $(LIBS)

testmainPathline: testmainPathline.o  lib$(LIBNAME).a
	$(C++) -o testmainPathline testmainPathline.o -L. -l$(LIBNAME) -lm 

testmainStreak: testmainStreak.o  lib$(LIBNAME).a
	$(C++) -o testmainStreak testmainStreak.o -L. -l$(LIBNAME) -lm

gldrawPathline: gldrawPathline.o  lib$(LIBNAME).a
	$(C++) -o gldrawPathline gldrawPathline.o -L. -l$(LIBNAME) $(LIBS)

gldrawPathline2: gldrawPathline2.o  lib$(LIBNAME).a
	$(C++) -o gldrawPathline2 gldrawPathline2.o -L. -l$(LIBNAME) $(LIBS)

gldrawPathline3: gldrawPathline3.o  lib$(LIBNAME).a
	$(C++) -o gldrawPathline3 gldrawPathline3.o -L. -l$(LIBNAME) $(LIBS)

gldrawPathline4: gldrawPathline4.o  lib$(LIBNAME).a
	$(C++) -o gldrawPathline4 gldrawPathline4.o -L. -l$(LIBNAME) $(LIBS)

gldrawFlash: gldrawFlash.o  lib$(LIBNAME).a
	$(C++) -o gldrawFlash gldrawFlash.o -L. -l$(LIBNAME) $(LIBS)

gldrawFlash2: gldrawFlash2.o  lib$(LIBNAME).a
	$(C++) -o gldrawFlash2 gldrawFlash2.o -L. -l$(LIBNAME) $(LIBS)

gldrawFlash3: gldrawFlash3.o  lib$(LIBNAME).a
	$(C++) -o gldrawFlash3 gldrawFlash3.o -L. -l$(LIBNAME) $(LIBS)

gldrawFlash4: gldrawFlash4.o  lib$(LIBNAME).a
	$(C++) -o gldrawFlash4 gldrawFlash4.o -L. -l$(LIBNAME) $(LIBS)

gldrawFlashData: gldrawFlashData.o  lib$(LIBNAME).a
	$(C++) -o gldrawFlashData gldrawFlashData.o -L. -l$(LIBNAME) $(LIBS)

gldrawFlashTime: gldrawFlashTime.o  lib$(LIBNAME).a
	$(C++) -o gldrawFlashTime gldrawFlashTime.o -L. -l$(LIBNAME) $(LIBS)

gldrawFlashPathline: gldrawFlashPathline.o  lib$(LIBNAME).a
	$(C++) -o gldrawFlashPathline gldrawFlashPathline.o -L. -l$(LIBNAME) $(LIBS)

gldrawStreak: gldrawStreak.o  lib$(LIBNAME).a
	$(C++) -o gldrawStreak gldrawStreak.o -L. -l$(LIBNAME) $(LIBS)

#gldrawStreak2: gldrawStreak2.o  lib$(LIBNAME).a
#	$(C++) -o gldrawStreak2 gldrawStreak2.o -L. -l$(LIBNAME) $(LIBS)

gldrawStreak3: gldrawStreak3.o  lib$(LIBNAME).a
	$(C++) -o gldrawStreak3 gldrawStreak3.o -L. -l$(LIBNAME) $(LIBS)



clean:
	rm -f *.o *.a

depend: Makefile $(SRCS)
	makedepend -fMakefile $(CCFLAGS) $(INCLUDE) $(SRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

