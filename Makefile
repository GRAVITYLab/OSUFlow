#SHELL = /bin/sh

ARCH = MAC_OSX
#ARCH = LINUX

LIBNAME = OSUFlow
RM = rm 
AR = ar cq

TOP = ..

ifeq ($(ARCH),MAC_OSX)
C++ = g++
CC  = gcc
CCFLAGS = -g -c -DMAC_OSX 
LIBS  = -framework GLUT -framework OpenGL 
endif

ifeq ($(ARCH),LINUX)
CC = /soft/apps/packages/mpich2-1.0.7rc1-gcc/bin/mpicc
C++ = /soft/apps/packages/mpich2-1.0.7rc1-gcc/bin/mpicxx
CCFLAGS = -g -c -DMPI -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
LIBS = -lm -lglut -lGL
endif

INCLUDE = -I.

OBJS =  Candidate.o  Grid.o  polynomials.o  TimeVaryingFieldLine.o \
	eigenvals.o  Interpolator.o  Rake.o	    Topology.o \
	eigenvecs.o  IsoSurf.o	     Solution.o     triangulator.o \
	Element.o    StreakLine.o    VectorMatrix.o Lattice.o \
	Field.o      PathLine.o      Streamline.o \
	FieldLine.o  Plot3DReader.o  TimeLine.o \
	OSUFlow.o calc_subvolume.o Lattice4D.o 

SRCS =  Candidate.C  Grid.C  polynomials.C  TimeVaryingFieldLine.C \
	eigenvals.C  Interpolator.C  Rake.C	    Topology.C \
	eigenvecs.C  IsoSurf.C	     Solution.C     triangulator.C \
	Element.C    StreakLine.C    VectorMatrix.C Lattice.C \
	Field.C      PathLine.C      Streamline.C \
	FieldLine.C  Plot3DReader.C  TimeLine.C \
	OSUFlow.C calc_subvoulme.C Lattice4D.C 

.SUFFIXES: .c .C

.c.o:
	$(CC) $(CCFLAGS) $(INCLUDE) -c $<

.C.o:
	$(C++) $(CCFLAGS) $(INCLUDE) -c $<

default: all

all: lib$(LIBNAME).a testmain testmain2 testmain3 testmain4 gldraw gldraw2 gldraw3 gldraw4 \
	testmainT testmainStreak gldrawT gldrawStreak  gldrawStk2  gldrawStk3 drawtest  mpitest

#all: lib$(LIBNAME).a 

lib$(LIBNAME).a : $(OBJS)
	$(RM) -f $@
#	/usr/lib/DCC/edg_prelink -v $(OBJS) 
	$(AR) $@ $(OBJS) 
#	cp $@ ../lib

testmain: testmain.o lib$(LIBNAME).a
	$(C++) -o testmain testmain.o -L. -l$(LIBNAME) -lm

testmain2: testmain2.o lib$(LIBNAME).a
	$(C++) -o testmain2 testmain2.o -L. -l$(LIBNAME) -lm

testmain3: testmain3.o lib$(LIBNAME).a
	$(C++) -o testmain3 testmain3.o -L. -l$(LIBNAME) -lm

testmain4: testmain4.o lib$(LIBNAME).a
	$(C++) -o testmain4 testmain4.o -L. -l$(LIBNAME) -lm

mpitest: MpiMain.o lib$(LIBNAME).a
	$(C++) -o mpitest MpiMain.o -L. -l$(LIBNAME) -lm


drawtest: MpiDraw.o lib$(LIBNAME).a
	$(C++) -o drawtest MpiDraw.o -L. -l$(LIBNAME) $(LIBS) 

gldraw: gldraw.o  lib$(LIBNAME).a
	$(C++) -o gldraw gldraw.o -L. -l$(LIBNAME) $(LIBS) 

gldraw2: gldraw2.o  lib$(LIBNAME).a
	$(C++) -o gldraw2 gldraw2.o -L. -l$(LIBNAME) $(LIBS) 

gldraw3: gldraw3.o  lib$(LIBNAME).a
	$(C++) -o gldraw3 gldraw3.o -L. -l$(LIBNAME) $(LIBS)

gldraw4: gldraw4.o  lib$(LIBNAME).a
	$(C++) -o gldraw4 gldraw4.o -L. -l$(LIBNAME) $(LIBS)

testmainT: testmainT.o  lib$(LIBNAME).a
	$(C++) -o testmainT testmainT.o -L. -l$(LIBNAME) -lm 

testmainStreak: testmainStreak.o  lib$(LIBNAME).a
	$(C++) -o testmainStreak testmainStreak.o -L. -l$(LIBNAME) -lm

gldrawT: gldrawT.o  lib$(LIBNAME).a
	$(C++) -o gldrawT gldrawT.o -L. -l$(LIBNAME) $(LIBS)

gldrawStreak: gldrawStreak.o  lib$(LIBNAME).a
	$(C++) -o gldrawStreak gldrawStreak.o -L. -l$(LIBNAME) $(LIBS)

gldrawStk2: gldrawStk2.o  lib$(LIBNAME).a
	$(C++) -o gldrawStk2 gldrawStk2.o -L. -l$(LIBNAME) $(LIBS)

gldrawStk3: gldrawStk3.o  lib$(LIBNAME).a
	$(C++) -o gldrawStk3 gldrawStk3.o -L. -l$(LIBNAME) $(LIBS)

clean:
	rm -f *.o *.a

depend: Makefile $(SRCS)
	makedepend -fMakefile $(CCFLAGS) $(INCLUDE) $(SRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

