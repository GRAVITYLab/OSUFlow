#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([BIL], [0.6.0], [kendall@eecs.utk.edu])
AC_CONFIG_SRCDIR([src/bil.c])
AC_CONFIG_HEADER([config.h])

# Automake
AM_INIT_AUTOMAKE

# Disabling pnetcdf
AC_ARG_ENABLE(pnetcdf,
    [  --disable-pnetcdf  disable parallel netcdf support],
   [disable_pnetcdf=true])
AM_CONDITIONAL(DISABLE_PNETCDF, test x$disable_pnetcdf = xtrue)

# Disabling timing
AC_ARG_ENABLE(timing,
    [AS_HELP_STRING([--disable-timing], [disable timing support])],
   [disable_timing=true])
AM_CONDITIONAL(DISABLE_TIMING, test x$disable_timing = xtrue)

AC_ARG_WITH(pnetcdf,
		[AS_HELP_STRING([--with-pnetcdf=/path/to/implementation],
                    [installation prefix for parallel netcdf implementation])],
		[with_pnetcdf="$withval/include"], [with_pnetcdf="no"])
AC_SUBST(with_pnetcdf)

AC_ARG_WITH(mpi,
		[AS_HELP_STRING([--with-mpi=/path/to/implementation],
                    [installation prefix for mpi implementation])],
		[with_mpi="$withval"])

# Allow the MPICC environment variable to override configuration
if test ! ${MPICC} ; then
  AC_PATH_PROG([MPICC], [mpicc], [no], [$with_mpi/bin$PATH_SEPARATOR$PATH])
  AC_SUBST(MPICC)

  if test ${MPICC} = "no" ; then
    AC_MSG_ERROR([cannot find mpicc: try using the --with-mpi option to specify the base path of your MPI installation. If MPI is not installed, it can be downloaded at http://www.mcs.anl.gov/research/projects/mpich2/]) 
  fi
fi

if test ${with_pnetcdf} != "no" ; then
  CFLAGS="$CFLAGS -I$with_pnetcdf "
  CPPFLAGS="$CPPFLAGS -I$with_pnetcdf "
fi

if test ! ${CC}; then
CC="gcc"
fi

saveCC="${CC}"
CC=${MPICC}
AC_CHECK_HEADERS([pnetcdf.h], [], [
if test "x$disable_pnetcdf" != xtrue ; then
AC_MSG_ERROR([cannot find pnetcdf.h: try using the --with-pnetcdf option if Parallel netCDF is installed. If Parallel netCDF is not installed, it can be downloaded at http://trac.mcs.anl.gov/projects/parallel-netcdf])
fi ], [#include "mpi.h"])
CC="$saveCC"

# Libtool
AC_PROG_LIBTOOL
AC_CONFIG_MACRO_DIR([m4])
AC_DISABLE_SHARED

# Checks for programs.
AM_PROG_CC_C_O
#AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([inttypes.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_INT64_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
