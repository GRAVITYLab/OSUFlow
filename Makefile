SHELL = /bin/sh

LIBNAME = OSUFlow
RM = rm 
AR = ar cq

TOP = ..
C++ = g++
CC  = gcc

CCFLAGS = 

INCLUDE = -I.

OBJS =  Candidate.o  Grid.o  polynomials.o  TimeVaryingFieldLine.o \
	eigenvals.o  Interpolator.o  Rake.o	    Topology.o \
	eigenvecs.o  IsoSurf.o	     Solution.o     triangulator.o \
	Element.o    StreakLine.o   VectorMatrix.o \
	Field.o      PathLine.o      Streamline.o \
	FieldLine.o  Plot3DReader.o  TimeLine.o \
	OSUFlow.o 

SRCS =  Candidate.C  Grid.C  polynomials.C  TimeVaryingFieldLine.C \
	eigenvals.C  Interpolator.C  Rake.C	    Topology.C \
	eigenvecs.C  IsoSurf.C	     Solution.C     triangulator.C \
	Element.C    StreakLine.C   VectorMatrix.C \
	Field.C      PathLine.C      Streamline.C \
	FieldLine.C  Plot3DReader.C  TimeLine.C \
	OSUFlow.C


.SUFFIXES: .c .C

.c.o:
	$(C++) $(CCFLAGS) $(INCLUDE) -c $<

.C.o:
	$(C++) $(CCFLAGS) $(INCLUDE) -c $<

default: all

all: lib$(LIBNAME).a testmain
#all: lib$(LIBNAME).a 

lib$(LIBNAME).a : $(OBJS)
	$(RM) -f $@
#	/usr/lib/DCC/edg_prelink -v $(OBJS) 
	$(AR) $@ $(OBJS) 
#	cp $@ ../lib

testmain: testmain.o lib$(LIBNAME).a
	$(C++) -o testmain testmain.o -L. -l$(LIBNAME) -lm

clean:
	rm -f *.o *.a

depend: Makefile $(SRCS)
	makedepend -fMakefile $(CCFLAGS) $(INCLUDE) $(SRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

