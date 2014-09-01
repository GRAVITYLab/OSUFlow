# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (gcb) 

# ADD-BY-LEETEN 08/23/2012-BEGIN
if( DEFINED WIN32 )
	message(WARNING "For Win32 environment, please manually specify the GLUT library to used (glut32,freeglut,openglut,etc)")
endif()
# ADD-BY-LEETEN 08/23/2012-END

# specify the path to the headers
include_directories(
		. 
	)

# grab all the header files
file(GLOB HEADER_FILES "*.h" )

# grab all the source codes
AUX_SOURCE_DIRECTORY ("." 	LIBGCB_SRC_FILES)

# build the project as a library
add_library(${PROJECT_NAME} ${LIBGCB_SRC_FILES} ${HEADER_FILES})

# ADD-BY-LEETEN 09/01/2011-BEGIN
set_target_properties(${PROJECT_NAME} PROPERTIES 
	DEBUG_OUTPUT_NAME 	"${PROJECT_NAME}_d"
	RELEASE_OUTPUT_NAME "${PROJECT_NAME}"
	)

install(
	FILES ${HEADER_FILES}
	DESTINATION include/osuflow/renderer/gcb)
	
install(
	TARGETS ${PROJECT_NAME}
	DESTINATION lib/osuflow)
# ADD-BY-LEETEN 09/01/2011-END

#
# $Log$
#

