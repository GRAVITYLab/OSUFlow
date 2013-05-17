# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (VTK_examples) 

find_package(OpenGL)
find_package(GLUT)


# specify the path to the headers
include_directories(
	${OSUFLOWVTK_DIR}
)

link_directories (
	${OSUFLOWVTK_DIR}
)

add_executable(StreamlineSimpleVTK StreamlineSimpleVTK.cpp)
target_link_libraries ( StreamlineSimpleVTK ${VTK_LIBRARIES} OSUFlow OSUFlowVTK )

if( ${WITH_GRAPHICS} )

	add_executable(glStreamlineSimpleVTK glStreamlineSimpleVTK.cpp)
	target_link_libraries ( glStreamlineSimpleVTK ${VTK_LIBRARIES} OSUFlow OSUFlowVTK GL GLU glut)

	# GCB related application
	include_directories ( ${RENDERER_DIR} ${GCB_DIR} )
	link_directories ( ${RENDERER_DIR} ${GCB_DIR} )

	add_executable(gcbStreamlineVTK gcbStreamlineVTK.cpp)
	target_link_libraries( gcbStreamlineVTK ${VTK_LIBRARIES} OSUFlow OSUFlowVTK gcb renderer GL GLU glut  )

	add_executable(StreamlineSimpleVTKPipeline StreamlineSimpleVTKPipeline.cpp)
	target_link_libraries( StreamlineSimpleVTKPipeline OSUFlowVTK ${VTK_LIBRARIES} OSUFlow GL GLU glut )

endif( ${WITH_GRAPHICS} )

