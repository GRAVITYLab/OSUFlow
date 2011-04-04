#----------------------------------------------------------------------------
#
# makefile system definitions
#
# Tom Peterka
# Argonne National Laboratory
# 9700 S. Cass Ave.
# Argonne, IL 60439
# tpeterka@mcs.anl.gov
#
# All rights reserved. May not be used, modified, or copied
# without permission
#
#----------------------------------------------------------------------------
#
# users: set your architecture, options, and paths in user_defs.h
# you should not need to touch this or other makefiles in the project
#
#----------------------------------------------------------------------------

# parallel netcdf, parallel hdf5, and bil don't make sense without MPI

ifneq ($(MPI), YES)
PNETCDF = NO
HDF5 = NO
BIL = NO
endif

# includes

INCLUDE =
ifeq ($(PNETCDF), YES)
INCLUDE += -I$(NETCDF_INC)
endif
ifeq ($(HDF5), YES)
INCLUDE += -I$(HDF_INC)
endif
ifeq ($(ZOLTAN), YES)
INCLUDE += -I$(ZOLTAN_INC)
endif
INCLUDE += $(MISC_INC) # no extra symbols such as -I are prepended for MISC_INC

# libraries

LIBS = -lm
ifeq ($(PNETCDF), YES)
LIBS += -L$(NETCDF_LIB)
endif
ifeq ($(ARCH), BGP)
ifeq ($(HDF5), YES)
LIBS += -L$(HDF_LIB)
endif
endif
ifneq ($(ARCH), BGP)
ifeq ($(HDF5), YES)
LIBS += -lz -L$(HDF_LIB)
endif
endif
ifeq ($(ZOLTAN), YES)
LIBS += -L$(ZOLTAN_LIB)
endif
ifeq ($(ARCH), MAC_OSX_OMPI)
ifeq ($(GRAPHICS), YES)
LIBS += -framework GLUT -framework OpenGL
endif
endif
ifeq ($(ARCH), MAC_OSX_MPICH)
ifeq ($(GRAPHICS), YES)
LIBS += -framework GLUT -framework OpenGL
endif
endif
ifeq ($(ARCH), LINUX)
ifeq ($(GRAPHICS), YES)
LIBS += -lglut -lGLU -lGL
endif
endif
ifeq ($(ARCH), LINUX_SERIAL)
ifeq ($(GRAPHICS), YES)
LIBS += -lglut -lGLU -lGL
endif
endif
LIBS += $(MISC_LIB) # no extra symbols such as -l are prepended for MISC_LIB

# compiler flags

CCFLAGS = -g -c
CCFLAGS += -D_OSUFLOW

# additional compiler flags and paths for each architecture

ifeq ($(ARCH), MAC_OSX_OMPI) # mac osx w/ open mpi
C++ = g++
CC  = gcc
CCFLAGS += -DMAC_OSX_OMPI
ifeq ($(MPI), YES)
CCFLAGS += -D_MPI
endif
endif

ifeq ($(ARCH), MAC_OSX_MPICH) # mac osx w/ mpich
C++ = mpicxx
CC = mpicc
CCFLAGS += -DMAC_OSX_MPICH
ifeq ($(MPI), YES)
CCFLAGS += -D_MPI
endif
endif

ifeq ($(ARCH), LINUX_SERIAL) # linux serial (used by Han-Wei)
C++ = g++
CC  = gcc
CCFLAGS += -DLINUX
LIBS =
endif

ifeq ($(ARCH), LINUX) # linux generic
C++ = mpicxx
CC = mpicc
ifeq ($(MPE), YES)
C++   = mpecxx -mpilog
endif
CCFLAGS += -DLINUX
ifeq ($(MPI), YES)
CCFLAGS += -D_MPI -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
endif
endif

ifeq ($(ARCH), BGP) # BG/P
C++ = mpixlcxx_r
CC = mpixlc_r
ifeq ($(MPE), YES)
C++ = /soft/apps/mpe/install_ibm/bin/mpecc -mpicc=mpixlcxx_r -mpilog
endif
CCFLAGS += -O3 -qarch=450d -qtune=450 -DBGP
ifeq ($(MPI), YES)
CCFLAGS += -D_MPI -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
endif
INCLUDE += -I/usr/local/include \
	-I/usr/X11R6/include \
	-I/bgsys/drivers/ppcfloor/arch/include \
	-I/soft/apps/hdf5-1.8.0/include \
	-I/soft/apps/parallel-netcdf-1.0.3-xl/include
endif

ifeq ($(ARCH), EUREKA) # eureka
INCLUDE = -I/usr/include
C++   = mpicxx
CC = mpicc
ifeq ($(MPE), YES)
C++   = mpecxx -mpilog
endif
CCFLAGS += -O3 -DEUREKA
ifeq ($(MPI), YES)
CCFLAGS += -D_MPI -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
endif
endif

# additional compiler flags

ifeq ($(PNETCDF), YES)
CCFLAGS += -DPNETCDF
endif
ifeq ($(HDF5), YES)
CCFLAGS += -DHDF5
endif
ifeq ($(GRAPHICS), YES)
CCFLAGS += -DGRAPHICS
endif
ifeq ($(MPE), YES)
CCFLAGS += -DMPE
endif
ifeq ($(BYTE_SWAP), YES)
CCFLAGS += -DBYTE_SWAP
endif
ifeq ($(DEBUG), YES)
CCFLAGS += -DDEBUG
endif
ifeq ($(WARNINGS), YES)
CCFLAGS += -Wall -Wextra
endif
ifeq ($(ZOLTAN), YES)
CCFLAGS += -DZOLTAN
endif
ifeq ($(BIL), YES)
CCFLAGS += -DUSE_BIL
endif
