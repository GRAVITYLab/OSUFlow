# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (OSUFlowVTK) 

set(WITH_PARAVIEW_PLUGINS	OFF	CACHE	BOOL	"Use Paraview Package")
set(OSUFLOWVTK_DIR      "${CMAKE_PROJECT_DIR}")


if (${WITH_MPI})
	find_package(MPI REQUIRED)
	include_directories(${MPI_CXX_INCLUDE_PATH})
endif (${WITH_MPI})

if (${WITH_GRAPHICS})
 find_package(OpenGL)
 find_package(GLUT)
endif (${WITH_GRAPHICS})



if (${WITH_PARAVIEW_PLUGINS})
	find_package(ParaView REQUIRED)
	include (${PARAVIEW_USE_FILE})
	# include (ParaViewBranding)
	include (ParaViewPlugins)
	# should enable argument -fPIC

	# Build VTK rendering modules when using Paraview package
	add_definitions(-DvtkRenderingCore_AUTOINIT=4\(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL\))
	add_definitions(-DvtkRenderingVolume_AUTOINIT=1\(vtkRenderingVolumeOpenGL\))

else()
	find_package(VTK REQUIRED)
	include (${VTK_USE_FILE})
endif()


# specify the path to the headers
include_directories(
	. 
	${OSUFLOW_DIR}/src
        ${OSUFLOWVTK_DIR}

)

# grab all the header files
file(GLOB HEADER_FILES "*.h" )

list(APPEND SRC_FILES
	VectorFieldVTK.cxx
	vtkOSUFlow.cxx
	)

# build the project as a library
add_library(${PROJECT_NAME} ${LIB_TYPE}  ${SRC_FILES} ${HEADER_FILES})

set_target_properties(${PROJECT_NAME} PROPERTIES 
        DEBUG_OUTPUT_NAME       "${PROJECT_NAME}_d"
        RELEASE_OUTPUT_NAME "${PROJECT_NAME}"
        )

install(
        FILES ${HEADER_FILES}
        DESTINATION include/VTK)

install(
        TARGETS ${PROJECT_NAME}
        DESTINATION lib)


#########################################################
# setup for dependent codes:
#########################################################

# setup link libraries
#message(STATUS "VTKLIBS: ${PARAVIEW_LIBRARIES}")
list (APPEND LINK_LIBS
	OSUFlowVTK
	OSUFlow
	${PARAVIEW_LIBRARIES}
	${VTK_LIBRARIES}
	)
if (WITH_GRAPHICS)
	include_directories ( 
		${OSUFLOW_DIR}/renderer 
		${OSUFLOW_DIR}/renderer/libgcb
	)
	list (APPEND LINK_LIBS
		renderer
		gcb
	)
        list (APPEND LINK_LIBS ${OPENGL_LIBRARY} ${GLU_LIBRARY} ${GLUT_LIBRARY})
endif ()

if (WITH_MPI)
	list (APPEND LINK_LIBS diy bil ${MPI_C_LIBRARIES})
endif ()
if (WITH_PNETCDF)
	list (APPEND LINK_LIBS  pnetcdf )
endif ()

	

# examples
add_subdirectory (examples)

# plugins
if (${WITH_PARAVIEW_PLUGINS})
	add_subdirectory (plugins)
endif()

# tests
if( ${WITH_GRAPHICS} )
	add_executable( compareStreamline compareStreamline.cpp )
        target_link_libraries ( compareStreamline ${LINK_LIBS})
endif ()
if ( ${WITH_MPI} AND ${WITH_PARAVIEW_PLUGINS})
	add_executable( testMPI testMPI.cpp )
	target_link_libraries ( testMPI ${LINK_LIBS} )
endif ()
