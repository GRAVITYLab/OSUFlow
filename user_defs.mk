#----------------------------------------------------------------------------
#
# makefile user definitions
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
# users: set your architecture, options, and paths in this file only
# you should not need to touch the other makefiles in the project
#
#----------------------------------------------------------------------------

# 1. Set your architecture here

#ARCH = MAC_OSX_OMPI
#ARCH = MAC_OSX_MPICH
ARCH = LINUX
#ARCH = LINUX_SERIAL
#ARCH = BGP
#ARCH = EUREKA

# 2. Set your options here

MPI = YES
MPE = NO
BIL = NO
PNETCDF = NO
HDF5 = NO
ZOLTAN = NO
GRAPHICS = YES
DEBUG = NO
BYTE_SWAP = NO
WARNINGS = NO

# 3. Set your paths here 
# keep a copy of them below for later reuse, commenting others' out
# please don't delete other people's settings

# Tom's linux
#HDF_INC = /homes/tpeterka/hdf5-install/include
#ZOLTAN_INC = /homes/tpeterka/Zoltan_v3.2/BuildDir/include
#NETCDF_INC = /homes/tpeterka/parallel-netcdf-1.0.2/src/lib
#MISC_INC =
#HDF_LIB = /homes/tpeterka/hdf5-install/lib
#ZOLTAN_LIB = /homes/tpeterka/Zoltan_v3.2/BuildDir/lib
#NETCDF_LIB = /homes/tpeterka/parallel-netcdf-1.0.2/src/lib
#MISC_LIB =

# Tom's mac
HDF_INC = /Users/tpeterka/software/hdf5-1.8.5/include
ZOLTAN_INC = /Users/tpeterka/software/zoltan-3.3/include
NETCDF_INC = /Users/tpeterka/software/parallel-netcdf-1.1.1/include
MISC_INC = 
HDF_LIB = /Users/tpeterka/software/hdf5-1.8.5/lib
ZOLTAN_LIB = /Users/tpeterka/software/zoltan-3.3/lib
NETCDF_LIB = /Users/tpeterka/software/parallel-netcdf-1.1.1/lib
MISC_LIB = 

# BG/P
#HDF_INC = /soft/apps/hdf5-1.8.0/include
#ZOLTAN_INC = 
#NETCDF_INC = /soft/apps/parallel-netcdf-1.0.3-xl/include
#MISC_INC = 
#HDF_LIB = /soft/apps/hdf5-1.8.0/lib
#ZOLTAN_LIB = 
#NETCDF_LIB = /soft/apps/parallel-netcdf-1.0.3-xl/lib
#MISC_LIB = 

# Boonth's linux
#HDF_INC = /home/nouanese/workspace/hdf5-1.8.4/build/include
#ZOLTAN_INC = /home/nouanese/workspace/Zoltan_v3.2/build/include
#NETCDF_INC = /home/nouanese/workspace/parallel-netcdf-1.1.1/build/include
#MISC_INC =
#HDF_LIB = /home/nouanese/workspace/hdf5-1.8.4/build/lib
#ZOLTAN_LIB = /home/nouanese/workspace/Zoltan_v3.2/build/lib
#NETCDF_LIB = /home/nouanese/workspace/parallel-netcdf-1.1.1/build/lib
#MISC_LIB = /home/nouanese/workspace/ParMETIS3_1/libparmetis.a /home/nouanese/workspace/ParMETIS3_1/libmetis.a 

# Han-Wei's mac
#HDF_INC = /Users/hwshen/Projects/hdf5-install/hdf5/include
#ZOLTAN_INC = 
#NETCDF_INC = 
#MISC_INC =
#HDF_LIB = /Users/hwshen/Projects/hdf5-install/hdf5/lib/libhdf5.a
#ZOLTAN_LIB = 
#NETCDF_LIB = 
#MISC_LIB =

# Han-Wei's linux serial
#HDF_INC = /home/hwshen/Libraries/hdf5-1.8.4-patch1/hdf5/include
#ZOLTAN_INC = 
#NETCDF_INC = 
#MISC_INC =
#HDF_LIB = /home/hwshen/Libraries/hdf5-1.8.4-patch1/hdf5/lib/libhdf5.a
#ZOLTAN_LIB = 
#NETCDF_LIB = 
#MISC_LIB =

# Wes' mac
#HDF_INC = 
#ZOLTAN_INC = 
#NETCDF_INC = 
#MISC_INC =
#HDF_LIB = 
#ZOLTAN_LIB = 
#NETCDF_LIB = 
#MISC_LIB =
