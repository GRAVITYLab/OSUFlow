# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (osuflow)

# # set paths to OSUFlow-related libraries
# set(OSUFLOW_DIR 	"${CMAKE_SOURCE_DIR}" 		CACHE	PATH "The path to OSUFlow.")
# set(RENDERER_DIR	"${OSUFLOW_DIR}/renderer" 	CACHE	PATH "The path to the renderer of OSUFlow.")
# set(GCB_DIR 		"${RENDERER_DIR}/libgcb" 	CACHE	PATH "The path to gcb (a GLUT-related mini window managemenet library).")

# Recurse into the "src" subdirectory. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (src)

add_subdirectory (renderer)
add_subdirectory (renderer/libgcb)
add_subdirectory (renderer/examples)

#
# $Log$
#

