# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (renderer) 

find_package(OpenGL)
find_package(GLUT)

# set paths to OSUFlow-related libraries
set(OSUFLOW_DIR 	"${CMAKE_SOURCE_DIR}" ) #		CACHE	PATH "The path to OSUFlow.")

if(NOT DEFINED WIN32 AND NOT DEFINED MAC_OSX_OMPI AND NOT DEFINED MAC_OSX_MPICH)
	add_definitions(
		-DLINUX
	)
endif()

# specify the path to the headers
include_directories(
		. 
		${OSUFLOW_DIR}/src
		${OPENGL_INCLUDE_PATH}
		${GLU_INCLUDE_PATH}
		${GLUT_INCLUDE_PATH}
	)

# grab all the header files
file(GLOB HEADER_FILES "*.h" )

list(APPEND SRC_FILES
	# ADD-BY-LEETEN 01/06/2012-BEGIN
	# Draw.C # Jimmy: this may not be needed when leeten_Draw.C is used			
	leeten_Draw.C 
	# ADD-BY-LEETEN 01/06/2012-END
	Renderer.cpp
	LineAnimatorInOpenGL.cpp       
	LineRenderer.cpp
	LineRendererInOpenGL.cpp       
	LineRenderer_ColorScheme.cpp
	TimeLineRendererInOpenGL.cpp   
	TubeRenderer.cpp
	TubeRendererInOpenGL.cpp
	)
	
# build the project as a library
add_library(${PROJECT_NAME} ${SRC_FILES} ${HEADER_FILES})

# ADD-BY-LEETEN 09/01/2011-BEGIN
set_target_properties(${PROJECT_NAME} PROPERTIES 
	DEBUG_OUTPUT_NAME 	"${PROJECT_NAME}_d"
	RELEASE_OUTPUT_NAME "${PROJECT_NAME}_r"
	)

install(
	FILES ${HEADER_FILES}
	DESTINATION include/renderer)
	
install(
	TARGETS ${PROJECT_NAME}
	DESTINATION lib)
# ADD-BY-LEETEN 09/01/2011-END

# subdirectories
add_subdirectory (libgcb)

if( ${WITH_EXAMPLES} ) # ADD-BY-LEETEN 12/16/2011
	add_subdirectory (examples)
endif ()

#
# $Log$
#

