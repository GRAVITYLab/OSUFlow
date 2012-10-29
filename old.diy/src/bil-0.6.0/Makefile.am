AUTOMAKE_OPTIONS = subdir-objects
lib_LIBRARIES = libbil.a
libbil_a_SOURCES = src/bil.c src/bil.h src/bil_misc.c src/bil_misc.h \
	src/bil_sched.c src/bil_sched.h src/bil_block.c src/bil_block.h \
	src/bil_pio.c src/bil_pio.h src/bil_timing.c src/bil_timing.h
include_HEADERS = src/bil.h
libbil_a_CFLAGS = $(AM_CFLAGS)
if DISABLE_PNETCDF
libbil_a_CFLAGS += -DDISABLE_PNETCDF
endif
if DISABLE_TIMING
libbil_a_CFLAGS += -DDISABLE_TIMING
endif
EXTRA_DIST = examples/makefile examples/README examples/test_bil_nc.c \
	examples/test_bil_raw.c examples/create_raw_files.c \
	examples/create_nc_files.c
CC=$(MPICC)
