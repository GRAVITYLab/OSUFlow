#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.61])
AC_INIT([diy], [0.1], [tpeterka@mcs.anl.gov])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

#---------------------------
#
# mpi
#
if test ! ${MPICC} ; then
  AC_PATH_PROG([MPICC], [mpicc], [no], [$with_mpi/bin$PATH_SEPARATOR$PATH])
  AC_PATH_PROG([MPICXX], [mpicxx], [no], [$with_mpi/bin$PATH_SEPARATOR$PATH])
  AC_SUBST(MPICC)
  AC_SUBST(MPICXX)

  if test ${MPICC} = "no" ; then
    AC_MSG_ERROR([cannot find mpicc: try using the --with-mpi option to specify the base path of your MPI installation. If MPI is not installed, it can be downloaded at http://www.mcs.anl.gov/research/projects/mpich2/]) 
  fi
fi

#---------------------------
#
# pnetcdf
#
AC_ARG_WITH(pnetcdf,
		[AS_HELP_STRING([--with-pnetcdf=/path/to/implementation],
                    [installation prefix for parallel netcdf implementation])],
		[with_pnetcdf="YES"], [with_pnetcdf="NO"])

if test ${with_pnetcdf} = "NO" ; then # path not given, if not automatic, disable it

  saveCC="${CC}"
  CC=${MPICC}
  AC_CHECK_HEADERS([pnetcdf.h], [with_pnetcdf="YES"])

  if test ${with_pnetcdf} = "NO" ; then
    AC_MSG_NOTICE([
    ----------------
    Cannot find parallel netcdf:  the package will be built without parallel netcdf support. Otherwise use --with-pnetcdf=/path/to/implementation to specify the base path (excluding lib or include subdirectories) of your parallel netcdf installation. Parallel netcdf can be downloaded at http://tracc.mcs.anl.gov/projects/parallel-netcdf
    ----------------]) 
  fi

  CC="$saveCC"

else # path given, check it

  saveCC="${CC}"
  CC=${MPICC}
  pnetcdf_inc="$withval/include"
  pnetcdf_lib="$withval/lib"
  CFLAGS="$CFLAGS -I$pnetcdf_inc "
  CPPFLAGS="$CPPFLAGS -I$pnetcdf_inc "
  AC_CHECK_HEADERS([pnetcdf.h], , [AC_MSG_ERROR([cannot find pnetcdf.h in the location $pnetcdf_inc. Double-check the base path provided in --with_pnetcdf argument? (Don't put /include or /lib at the end)])])
  AC_SUBST(pnetcdf_inc)
  AC_SUBST(pnetcdf_lib)

  CC="$saveCC"

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

if test ${with_hdf5} = "NO" ; then # path not given, if not automatic, disable it

  saveCC="${CC}"
  CC=${MPICC}
  AC_CHECK_HEADERS([hdf5.h], [with_hdf5="YES"])

  if test ${with_hdf5} = "NO" ; then
    AC_MSG_NOTICE([
    ----------------
    Cannot find parallel hdf5:  the package will be built without parallel hdf5 support. Otherwise use --with-hdf5=/path/to/implementation to specify the base path (excluding lib or include subdirectories) of your parallel hdf5 installation. Hdf5 can be downloaded at http://www.hdfgroup.org/HDF5
    ----------------]) 
  fi

  CC="$saveCC"

else # path given, check it

  saveCC="${CC}"
  CC=${MPICC}
  hdf5_inc="$withval/include"
  hdf5_lib="$withval/lib"
  CFLAGS="$CFLAGS -I$hdf5_inc "
  CPPFLAGS="$CPPFLAGS -I$hdf5_inc "
  AC_CHECK_HEADERS([hdf5.h], , [AC_MSG_ERROR([cannot find hdf5.h in the location $hdf5_inc. Double-check the base path provided in the --with_hdf5 argument? (Don't put /include or /lib at the end)])])
  AC_SUBST(hdf5_inc)
  AC_SUBST(hdf5_lib)

  CC="$saveCC"

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

if test ${with_zoltan} = "NO" ; then # path not given, if not automatic, disable it

  saveCC="${CC}"
  CC=${MPICC}
  AC_CHECK_HEADERS([zoltan.h], [with_zoltan="YES"])

  if test ${with_zoltan} = "NO" ; then
    AC_MSG_NOTICE([
    ----------------
    Cannot find zoltan:  the package will be built without zoltan support. Otherwise use --with-zoltan=/path/to/implementation to specify the base path (excluding lib or include subdirectories) of your zoltan installation. Zoltan can be downloaded at http://www.cs.sandia.gov/Zoltan
    ----------------]) 
  fi

  CC="$saveCC"

else # path given, check it

  saveCC="${CC}"
  CC=${MPICC}
  zoltan_inc="$withval/include"
  zoltan_lib="$withval/lib"
  CFLAGS="$CFLAGS -I$zoltan_inc "
  CPPFLAGS="$CPPFLAGS -I$zoltan_inc "
  AC_CHECK_HEADERS([zoltan.h], , [AC_MSG_ERROR([cannot find zoltan.h in the location $zoltan_inc. Double-check the base path provided in the --with_zoltan argument? (Don't put /include or /lib at the end)])])
  AC_SUBST(zoltan_inc)
  AC_SUBST(zoltan_lib)

  CC="$saveCC"

fi

echo "zoltan enabled = $with_zoltan"
AC_SUBST(with_zoltan)

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
# Checks for libraries.
#
AC_CHECK_LIB([m], [pow])

# Checks for header files.
AC_CHECK_HEADERS([stddef.h stdint.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# configures
AC_CONFIG_FILES([config_defs
                 makefile
                 src/makefile])

if test ${enable_bil} = "YES"; then
AC_CONFIG_SUBDIRS([src/bil-0.6.0])
fi

AC_OUTPUT
