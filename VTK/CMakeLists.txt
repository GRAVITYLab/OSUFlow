# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (OSUFlowVTK) 

find_package(VTK 6 REQUIRED)
include (${VTK_USE_FILE})

set(OSUFLOWVTK_DIR	"${CMAKE_PROJECT_DIR}")

# specify the path to the headers
include_directories(
	. 
	${OSUFLOW_DIR}/src
        ${OSUFLOWVTK_DIR}

)

# grab all the header files
file(GLOB HEADER_FILES "*.h" )

list(APPEND SRC_FILES
	vtkOSUFlow.cxx
	VectorFieldVTK.cxx
	OSUFlowVTK.cxx
	)
	
# build the project as a library
add_library(${PROJECT_NAME} ${SRC_FILES} ${HEADER_FILES})

set_target_properties(${PROJECT_NAME} PROPERTIES 
	DEBUG_OUTPUT_NAME 	"${PROJECT_NAME}_d"
	RELEASE_OUTPUT_NAME "${PROJECT_NAME}"
	)

install(
	FILES ${HEADER_FILES}
	DESTINATION include/VTK)
	
install(
	TARGETS ${PROJECT_NAME}
	DESTINATION lib)


# examples
add_subdirectory (examples)

# tests
if( ${WITH_GRAPHICS} )

	add_executable( compareStreamline compareStreamline.cpp )
        target_link_libraries ( compareStreamline OSUFlowVTK ${VTK_LIBRARIES} OSUFlow GL GLU glut)

endif ()
