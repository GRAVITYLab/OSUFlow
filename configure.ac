#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.61])
AC_INIT([osuflow], [1.0], [tpeterka@mcs.anl.gov])
AC_CONFIG_SRCDIR([src/ParFlow.C])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# set user options

#---------------------------
#
# architecture
#
AC_ARG_WITH(arch,
		[AS_HELP_STRING([--with-arch=<architecture>],
[Architecture is a mandatory setting. Valid architectures are (case sensitive): MAC_OSX_OMPI, MAC_OSX_MPICH, LINUX, LINUX_SERIAL, BGP, EUREKA])], 
[if test ${withval} != "MAC_OSX_OMPI" && ${withval} != "MAC_OSX_MPICH" && ${withval} != "LINUX" && ${withval} != "LINUX_SERIAL" && ${withval} != "BGP" && ${withval} != "EUREKA" ; then
  AC_MSG_ERROR([Please specify an architectures using --with-arch=<architecture>. Valid architectures are (case sensitive): 
MAC_OSX_OMPI
MAC_OSX_MPICH
LINUX
LINUX_SERIAL
BGP
EUREKA])
fi],
[AC_MSG_ERROR([Please specify an architecture using --with-arch=<architecture>. Valid architectures are (case sensitive): 
MAC_OSX_OMPI
MAC_OSX_MPICH
LINUX
LINUX_SERIAL
BGP
EUREKA])])
echo "architecture = $withval"
AC_SUBST(arch, $withval)

#---------------------------
#
# mpe
#
AC_ARG_ENABLE(mpe, [AS_HELP_STRING([--enable-mpe],
[whether to enable mpe profiling for later viewing in jumpshot])],
[if test ${enableval} = "yes" ; then
  enable_mpe="YES"
else
  enable_mpe="NO"
fi],
[enable_mpe="NO"])
echo "mpe enabled = $enable_mpe"
AC_SUBST(enable_mpe)

#---------------------------
#
# bil
#
AC_ARG_ENABLE(bil, [AS_HELP_STRING([--disable-bil],
[whether to disable bil (block I/O library) from performing I/O operations (enabled by default)])],
[if test ${enableval} = "yes" ; then
  enable_bil="YES"
else
  enable_bil="NO"
fi],
[enable_bil="YES"])
echo "bil enabled = $enable_bil"
AC_SUBST(enable_bil)

#---------------------------
#
# graphics
#
if test $with_arch = "BGP" ; then
  enable_graphics="NO"
else
  AC_ARG_ENABLE(graphics, [AS_HELP_STRING([--disable-graphics],
    [whether to disable graphics (enabled by default)])],
    [if test ${enableval} = "yes" ; then
    enable_graphics="YES"
    else
    enable_graphics="NO"
    fi],
    [enable_graphics="YES"])
fi
echo "graphics enabled = $enable_graphics"
AC_SUBST(enable_graphics)

#---------------------------
#
# debug
#
AC_ARG_ENABLE(debug, [AS_HELP_STRING([--enable-debug],
[whether to enable debug mode (numerous verbose printfs)])],
[if test ${enableval} = "yes" ; then
  enable_debug="YES"
else
  enable_debug="NO"
fi],
[enable_debug="NO"])
echo "debug enabled = $enable_debug"
AC_SUBST(enable_debug)

#---------------------------
#
# byte swap
#
AC_ARG_ENABLE(byteswap, [AS_HELP_STRING([--enable-byteswap],
[whether to enable byte swapping (when running on a machine with opposite endianness than the dataset)])],
[if test ${enableval} = "yes" ; then
  enable_byteswap="YES"
else
  enable_byteswap="NO"
fi],
[enable_byteswap="NO"])
echo "byteswap enabled = $enable_byteswap"
AC_SUBST(enable_byteswap)

#---------------------------
#
# warnings
#
AC_ARG_ENABLE(warnings, [AS_HELP_STRING([--enable-warnings],
[whether to enable all warnings during build])],
[if test ${enableval} = "yes" ; then
  enable_warnings="YES"
else
  enable_warnings="NO"
fi],
[enable_warnings="NO"])
echo "warnings enabled = $enable_warnings"
AC_SUBST(enable_warnings)

#---------------------------
#
# mpi enabled
#
AC_ARG_ENABLE(mpi, [AS_HELP_STRING([--disable-mpi],
[whether to disable mpi (build serial mode only)])],
[if test ${enableval} = "no" ; then
  with_mpi="NO"
else
  with_mpi="YES"
fi],
[with_mpi="YES"])
echo "mpi enabled = $with_mpi"

#---------------------------
#
# mpi
#
if test $with_mpi = "YES" ; then

if test ! ${MPICC} ; then
  AC_PATH_PROG([MPICC], [mpicc], [no], [$with_mpi/bin$PATH_SEPARATOR$PATH])
  AC_PATH_PROG([MPICXX], [mpicxx], [no], [$with_mpi/bin$PATH_SEPARATOR$PATH])
  AC_SUBST(MPICC)
  AC_SUBST(MPICXX)

  if test ${MPICC} = "no" ; then
    AC_MSG_ERROR([cannot find mpicc: try using the --with-mpi option to specify the base path of your MPI installation. If MPI is not installed, it can be downloaded at http://www.mcs.anl.gov/research/projects/mpich2/]) 
  fi
fi

fi

AC_SUBST(with_mpi)

#---------------------------
#
# pnetcdf
#
AC_ARG_WITH(pnetcdf,
		[AS_HELP_STRING([--with-pnetcdf=/path/to/implementation],
                    [installation prefix for parallel netcdf implementation])],
		[with_pnetcdf="YES"], [with_pnetcdf="NO"])
if test ${with_mpi} = "YES" ; then
if test ${with_pnetcdf} = "NO" ; then
  AC_CHECK_HEADERS([pnetcdf.h], [with_pnetcdf="YES"])
else
  pnetcdf_inc="$withval/include"
  pnetcdf_lib="$withval/lib"
  CFLAGS="$CFLAGS -I$pnetcdf_inc "
  CPPFLAGS="$CPPFLAGS -I$pnetcdf_inc "
  saveCC="${CC}"
  CC=${MPICC}
  AC_CHECK_HEADERS([pnetcdf.h], [], [AC_MSG_ERROR([cannot find pnetcdf.h in the location $pnetcdf_inc. Double-check the base path provided in --with-pnetcdf argument? (Don't put /include or /lib at the end)])], [#include "mpi.h"])
  CC="$saveCC"
  AC_SUBST(pnetcdf_inc)
  AC_SUBST(pnetcdf_lib)
fi
if test ${with_pnetcdf} = "NO" ; then
  AC_MSG_NOTICE([
----------------
Cannot find parallel netcdf:  the package will be built without parallel netcdf support. Otherwise use --with-pnetcdf=/path/to/implementation to specify the base path (excluding lib or include subdirectories) of your parallel netcdf installation. Parallel netcdf can be downloaded at http://tracc.mcs.anl.gov/projects/parallel-netcdf
----------------]) 
fi
fi
echo "pnetcdf enabled = $with_pnetcdf"
AC_SUBST(with_pnetcdf)

#---------------------------
#
# hdf5
#
AC_ARG_WITH(hdf5,
		[AS_HELP_STRING([--with-hdf5=/path/to/implementation],
                    [installation prefix for parallel netcdf implementation])],
		[with_hdf5="YES"], [with_hdf5="NO"])
if test ${with_mpi} = "YES" ; then
if test ${with_hdf5} = "NO" ; then
  AC_CHECK_HEADERS([hdf5.h], [with_hdf5="YES"])
else
  hdf5_inc="$withval/include"
  hdf5_lib="$withval/lib"
  CFLAGS="$CFLAGS -I$hdf5_inc "
  CPPFLAGS="$CPPFLAGS -I$hdf5_inc "
  saveCC="${CC}"
  CC=${MPICC}
  AC_CHECK_HEADERS([hdf5.h], , [AC_MSG_ERROR([cannot find hdf5.h in the location $hdf5_inc. Double-check the base path provided in the --with_hdf5 argument? (Don't put /include or /lib at the end)])], [#include "mpi.h"])
  AC_SUBST(hdf5_inc)
  AC_SUBST(hdf5_lib)
  CC="$saveCC"
fi
if test ${with_hdf5} = "NO" ; then
  AC_MSG_NOTICE([
----------------
Cannot find parallel hdf5:  the package will be built without parallel hdf5 support. Otherwise use --with-hdf5=/path/to/implementation to specify the base path (excluding lib or include subdirectories) of your parallel hdf5 installation. Hdf5 can be downloaded at http://www.hdfgroup.org/HDF5
----------------]) 
fi
fi
echo "hdf5 enabled = $with_hdf5"
AC_SUBST(with_hdf5)

#---------------------------
#
# zoltan
#
AC_ARG_WITH(zoltan,
		[AS_HELP_STRING([--with-zoltan=/path/to/implementation],
                    [installation prefix for parallel netcdf implementation])],
		[with_zoltan="YES"], [with_zoltan="NO"])
if test ${with_mpi} = "YES" ; then
if test ${with_zoltan} = "NO" ; then
  AC_CHECK_HEADERS([zoltan.h], [with_zoltan="YES"])
else
  zoltan_inc="$withval/include"
  zoltan_lib="$withval/lib"
  CFLAGS="$CFLAGS -I$zoltan_inc "
  CPPFLAGS="$CPPFLAGS -I$zoltan_inc "
  AC_CHECK_HEADERS([zoltan.h], , [AC_MSG_ERROR([cannot find zoltan.h in the location $zoltan_inc. Double-check the base path provided in the --with_zoltan argument? (Don't put /include or /lib at the end)])])
  AC_SUBST(zoltan_inc)
  AC_SUBST(zoltan_lib)
fi
if test ${with_zoltan} = "NO" ; then
  AC_MSG_NOTICE([
----------------
Cannot find zoltan:  the package will be built without zoltan support. Otherwise use --with-zoltan=/path/to/implementation to specify the base path (excluding lib or include subdirectories) of your zoltan installation. Zoltan can be downloaded at http://www.cs.sandia.gov/Zoltan
----------------]) 
fi
fi
echo "zoltan enabled = $with_zoltan"
AC_SUBST(with_zoltan)

#---------------------------
#
# build shared object in addition to static lib
#
AC_ARG_ENABLE(so, [AS_HELP_STRING([--enable-so],
[whether to generate dynamically-linked shared object in addition to statically-linked library (disabled by default)])],
[if test ${enableval} = "yes" ; then
  enable_so="YES"
else
  enable_so="NO"
fi],
[enable_so="NO"])

echo "shared object enabled = $enable_so"
AC_SUBST(enable_so)

#---------------------------

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h float.h limits.h memory.h stdlib.h string.h strings.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([atexit floor memmove memset pow sqrt strchr strdup strrchr])

AC_CONFIG_FILES([AMR/Parallel/Makefile
                 AMR/Serial/Makefile
                 Makefile
		 config_defs
                 Regular/Parallel/Makefile
                 Regular/Serial/Makefile
                 src/Makefile
                 renderer/Makefile
                 renderer/examples/Makefile
                 renderer/libgcb/Makefile])

if test ${with_mpi} = "YES"; then
AC_CONFIG_SUBDIRS([diy])
fi

AC_OUTPUT
