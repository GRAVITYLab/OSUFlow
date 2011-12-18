# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (osuflow)

# # set paths to OSUFlow-related libraries
# set(OSUFLOW_DIR 	"${CMAKE_SOURCE_DIR}" 		CACHE	PATH "The path to OSUFlow.")
# set(RENDERER_DIR	"${OSUFLOW_DIR}/renderer" 	CACHE	PATH "The path to the renderer of OSUFlow.")
# set(GCB_DIR 		"${RENDERER_DIR}/libgcb" 	CACHE	PATH "The path to gcb (a GLUT-related mini window managemenet library).")

# ADD-BY-LEETEN 08/17/2011-BEGIN
set(WITH_MPI FALSE CACHE BOOL "Specify whether MPI is used")
set(DIY_DIR "/usr" CACHE PATH "The installed path of DIY")

set(WITH_BIL	FALSE CACHE BOOL "Specify whether BIL is used.")	# ADD-BY-LEETEN 12/17/2011
set(BIL_DIR "/usr" CACHE PATH "The installed path of BIL")

set(WITH_PNETCDF FALSE CACHE BOOL "Specify whether PNETCDF is used")
set(PNETCDF_PATH "/usr" CACHE PATH "The installed path of PNETCDF")

if(NOT DEFINED WIN32 AND NOT DEFINED MAC_OSX_OMPI AND NOT DEFINED MAC_OSX_MPICH)
	add_definitions(
		-DLINUX
	)
endif()

# ADD-BY-LEETEN 12/16/2011-BEGIN
set(WITH_UNISTD FALSE CACHE BOOL "Specify whether unistd.h is included.")
if( ${WITH_UNISTD} )
	add_definitions(
		-DWITH_UNISTD
	)
endif()
# ADD-BY-LEETEN 12/16/2011-END

add_definitions(
	-D_OSUFLOW
	)

set(OSUFLOW_DIR 	"${CMAKE_SOURCE_DIR}")
# MOD-BY-LEETEN 09/01/2011-FROM:
	# set(OSUFLOW_INC		"${OSUFLOW_DIR}/src")	
	# set(OSUFLOW_LIB		"${OSUFLOW_DIR}/src")	
# MOD-BY-LEETEN 09/01/2011-TO:
set(OSUFLOW_SRC		"${OSUFLOW_DIR}/src")	
set(OSUFLOW_LIB		"${OSUFLOW_DIR}/lib")	
set(OSUFLOW_INC		"${OSUFLOW_DIR}/include")	
# MOD-BY-LEETEN 09/01/2011-END
set(RENDERER_DIR	"${OSUFLOW_DIR}/renderer" 	)
set(GCB_DIR 		"${RENDERER_DIR}/libgcb" 	)

include_directories(
	# ADD-BY-LEETEN 09/01/2011-BEGIN
	${OSUFLOW_SRC}
	# ADD-BY-LEETEN 09/01/2011-END
	${OSUFLOW_INC}
)

link_directories(
	${OSUFLOW_LIB}
)

if( WITH_MPI )
	find_package(MPI REQUIRED)
	
	set(DIY_INC "${DIY_DIR}/include")
	set(DIY_LIB "${DIY_DIR}/lib")
	
	set(BIL_INC "${BIL_DIR}/include")
	set(BIL_LIB "${BIL_DIR}/lib")

	add_definitions(	
		-D_MPI
	)
	include_directories(
		${MPI_INCLUDE_PATH}
		${DIY_INC}
		${BIL_INC}
	)
	link_directories(
		${DIY_LIB}
		${BIL_LIB}
	)

	if( WITH_PNETCDF )
		set(PNETCDF_INC "${PNETCDF_PATH}/include") 
		set(PNETCDF_LIB "${PNETCDF_PATH}/lib") 
	
		# MOD-BY-LEETEN 12/16/2011-FROM:
		  # add_definitions(-IPNETCDF)
		# TO:
		add_definitions(-DPNETCDF)
		# MOD-BY-LEETEN 12/16/2011-END
		include_directories(
			${PNETCDF_INC}
		)
		link_directories(
			${PNETCDF_LIB}
		)
	endif( )
endif()
# ADD-BY-LEETEN 08/17/2011-END

# Recurse into the "src" subdirectory. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (src)

# ADD-BY-LEETEN 12/16/2011-BEGIN
set(WITH_EXAMPLES  TRUE CACHE BOOL "Specify whether sample programs should be built as well.") # MOD-BY-LEETEN 12/16/2011-FROM: set(WITH_EXAMPLE  TRUE CACHE BOOL "Specify whether sample programs should be built as well.")
set(WITH_GRAPHICS TRUE CACHE BOOL "Specify whether Graphics-related libraries should be built.")
if( ${WITH_GRAPHICS} )
# ADD-BY-LEETEN 12/16/2011-END

add_subdirectory (renderer)
add_subdirectory (renderer/libgcb)
		 if( ${WITH_EXAMPLES} ) # ADD-BY-LEETEN 12/16/2011
		     add_subdirectory (renderer/examples)

# ADD-BY-LEETEN 12/16/2011-BEGIN
		 endif() 
endif() 
# ADD-BY-LEETEN 12/16/2011-END

# ADD-BY-LEETEN 08/17/2011-BEGIN
if( ${WITH_MPI} ) # MOD-BY-LEETEN 12/16/2011-FROM: if( WITH_MPI )
    if( ${WITH_EXAMPLES} ) # ADD-BY-LEETEN 12/16/2011
	add_subdirectory (Regular/Parallel)
    endif() # ADD-BY-LEETEN 12/16/2011
endif()
# ADD-BY-LEETEN 08/17/2011-END

# ADD-BY-LEETEN 12/17/2011-BEGIN
set(ARCH	 "LINUX" CACHE STRING	"Specify the architecture.")
set(WITH_MPE FALSE CACHE BOOL "Specify whether MPE is used.")

set(WITH_ZOLTAN FALSE CACHE BOOL "Specify whether Zoltabn is used.")
set(ZOLTAN_DIR 	"/usr" CACHE PATH "The installed path of Zoltan.")

set(WITH_HDF5 	FALSE CACHE BOOL "Specify whether HDF5 is used.")
set(HDF5_DIR 	"/usr" CACHE PATH "The installed path of HDF5.")

set(IS_DEBUG 	FALSE CACHE BOOL "Specify whether DEBUG mode is on.")
set(IS_BYTE_SWAP 	FALSE CACHE BOOL "Specify whether the byte should be swapped.")
set(IS_WARNING 		FALSE CACHE BOOL "Specify whether the warning message should be dumpped.")

# write a macro to create new variables whose name begins with "FLAG_"
macro(BECOME_FLAG	VAR)
	set(FLAG_${VAR}	NO)
	if( ${VAR} )
		set(FLAG_${VAR}	YES)
	endif()
endmacro(BECOME_FLAG) 

BECOME_FLAG(WITH_MPI)
BECOME_FLAG(WITH_MPE)
BECOME_FLAG(WITH_BIL)
BECOME_FLAG(WITH_PNETCDF)
BECOME_FLAG(WITH_HDF5)
BECOME_FLAG(WITH_GRAPHICS)
BECOME_FLAG(WITH_ZOLTAN)
BECOME_FLAG(IS_DEBUG)
BECOME_FLAG(IS_BYTE_SWAP)
BECOME_FLAG(IS_WARNING)

configure_file (
  "${PROJECT_SOURCE_DIR}/user_defs.mk.in"
  "${PROJECT_SOURCE_DIR}/user_defs.mk"
  )
# ADD-BY-LEETEN 12/17/2011-END

#
# $Log$
#

