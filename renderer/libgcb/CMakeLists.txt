# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (gcb) 

find_package(OpenGL)
find_package(GLUT)

if(NOT DEFINED WIN32 AND NOT DEFINED MAC_OSX_OMPI AND NOT DEFINED MAC_OSX_MPICH)
	add_definitions(
		-DLINUX
	)
endif()

# specify the path to the headers
include_directories(
		. 
	)

# grab all the header files
file(GLOB HEADER_FILES "*.h" )

# grab all the source codes
AUX_SOURCE_DIRECTORY ("." 	SRC_FILES)

# build the project as a library
add_library(${PROJECT_NAME} ${SRC_FILES} ${HEADER_FILES})

#
# $Log$
#

