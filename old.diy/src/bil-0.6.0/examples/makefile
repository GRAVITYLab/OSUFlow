############################################################################
 #   Copyright (C) 2010  Wes Kendall                                       #
 #   kendall@eecs.utk.edu                                                  #
 #                                                                         #
 #   This program is free software; you can redistribute it and/or modify  #
 #   it under the terms of the GNU Lesser General Public License as        #
 #   published by the Free Software Foundation; either version 3 of the    #
 #   License, or (at your option) any later version.                       #
 #                                                                         #
 #   This program is distributed in the hope that it will be useful,       #
 #   but WITHOUT ANY WARRANTY; without even the implied warranty of        #
 #   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
 #   GNU Lesser General Public License for more details.                   #
 #                                                                         #
 #   You should have received a copy of the GNU Lesser General Public      #
 #   License along with this program; if not, write to the                 #
 #   Free Software Foundation, Inc.,                                       #
 #   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             #
 ############################################################################

# makefile
# Wes Kendall
# 05/29/2010
# simple makefile for compiling and running the examples. NOTE 
# that if you compiled BIL without netCDF support, you will 
# have to set DISABLE_NETCDF=1
####

# set to 1 if you disabled netcdf during configuration
DISABLE_NETCDF=1

INCLUDE=-I../src
CC=mpicc
CFLAGS=-Wall -O3

ifeq ($(DISABLE_NETCDF), 1)
	EXECS=test_bil_raw create_raw_files
	LIB=-L../ -lbil
else
	EXECS=test_bil_raw create_raw_files test_bil_nc create_nc_files
	LIB=-lnetcdf -L../ -lbil -lpnetcdf -lcurl
endif

.SUFFIXES: .c
.c.o:
	$(CC) $(CFLAGS) $(INCLUDE) -c $<

all: $(EXECS)

create_nc_files: create_nc_files.o
	$(CC) -o $@ create_nc_files.o $(LIB)

test_bil_nc: test_bil_nc.o create_nc_files
	$(CC) -o $@ test_bil_nc.o $(LIB)

create_raw_files: create_raw_files.o
	$(CC) -o $@ create_raw_files.o

test_bil_raw: test_bil_raw.o create_raw_files
	$(CC) -o $@ test_bil_raw.o $(LIB)

# arguments for running the example testing programs for BIL
TEST_FILE_DIR=./ # directory where testing files are stored
TEST_FILE_DIM_SIZE=16 # dimension size of each file
TEST_NUM_FILES=4 # number of files to generate for testing
TEST_NUM_PROCS=3 # number of procs to use for testing
ARGS=$(TEST_FILE_DIM_SIZE) $(TEST_NUM_FILES) $(TEST_FILE_DIR)

# runs the example raw reader. it first creates raw files and then
# runs the program with mpirun. for submitting job to a system, mpirun
# will have to be replaced with the appropriated queueing command. to
# save the test files to a different directory, update the TEST_FILE_DIR
# variable
run_test_bil_raw: test_bil_raw
	./create_raw_files $(ARGS)
	mpirun -n $(TEST_NUM_PROCS) ./test_bil_raw $(ARGS)

# runs the example netcdf reader. it first creates netcdf files and then
# runs the program with mpirun. for submitting job to a system, mpirun
# will have to be replaced with the appropriated queueing command. to
# save the test files to a different directory, update the TEST_FILE_DIR
# variable
run_test_bil_nc: test_bil_nc
	./create_nc_files $(ARGS)
	mpirun -n $(TEST_NUM_PROCS) ./test_bil_nc $(ARGS)

clean:
	rm -f *.o $(EXECS)
