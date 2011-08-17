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
set(BIL_DIR "/usr" CACHE PATH "The installed path of BIL")

set(WITH_PNETCDF FALSE CACHE BOOL "Specify whether PNETCDF is used")
set(PNETCDF_PATH "/usr" CACHE PATH "The installed path of PNETCDF")

if(NOT DEFINED WIN32 AND NOT DEFINED MAC_OSX_OMPI AND NOT DEFINED MAC_OSX_MPICH)
	add_definitions(
		-DLINUX
	)
endif()

add_definitions(
	-D_OSUFLOW
	)

set(OSUFLOW_DIR 	"${CMAKE_SOURCE_DIR}")
set(OSUFLOW_INC		"${OSUFLOW_DIR}/src")	
set(OSUFLOW_LIB		"${OSUFLOW_DIR}/src")	
set(RENDERER_DIR	"${OSUFLOW_DIR}/renderer" 	)
set(GCB_DIR 		"${RENDERER_DIR}/libgcb" 	)

include_directories(
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
	
		add_definitions(-IPNETCDF)
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

add_subdirectory (renderer)
add_subdirectory (renderer/libgcb)
add_subdirectory (renderer/examples)

# ADD-BY-LEETEN 08/17/2011-BEGIN
if( WITH_MPI )
	add_subdirectory (Regular/Parallel)
endif()
# ADD-BY-LEETEN 08/17/2011-END

#
# $Log$
#

