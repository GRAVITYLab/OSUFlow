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


# all executables
add_executable(StreamlineSimpleVTK StreamlineSimpleVTK.cpp)
target_link_libraries ( StreamlineSimpleVTK ${LINK_LIBS})

add_executable(testVTS testVTS.cpp)
target_link_libraries ( testVTS ${VTK_LIBRARIES} )

add_executable(vec2vts vec2vts.cpp )
target_link_libraries ( vec2vts ${VTK_LIBRARIES} )

add_executable(vec2vti vec2vti.cpp )
target_link_libraries ( vec2vti ${VTK_LIBRARIES} )


if( ${WITH_GRAPHICS} )

	add_executable(glStreamlineSimpleVTK glStreamlineSimpleVTK.cpp)
	target_link_libraries ( glStreamlineSimpleVTK ${LINK_LIBS})

	# GCB related application
	include_directories ( ${RENDERER_DIR} ${GCB_DIR} )
	link_directories ( ${RENDERER_DIR} ${GCB_DIR} )

	add_executable(gcbStreamlineVTK gcbStreamlineVTK.cpp)
	target_link_libraries( gcbStreamlineVTK  gcb renderer ${LINK_LIBS}  )

	add_executable(StreamlineSimpleVTKPipeline StreamlineSimpleVTKPipeline.cpp)
	target_link_libraries( StreamlineSimpleVTKPipeline  ${LINK_LIBS} )

	add_executable(StreamlineRenderVec StreamlineRenderVec.cpp)
	target_link_libraries( StreamlineRenderVec  ${LINK_LIBS} )

	add_executable(StreamlineRenderPlot3D StreamlineRenderPlot3D.cpp)
	target_link_libraries( StreamlineRenderPlot3D  ${LINK_LIBS})

endif( ${WITH_GRAPHICS} )

### MPI Examples ###
if ( WITH_MPI )
	add_executable(MpiDrawVTK MpiDrawVTK.C)
	target_link_libraries ( MpiDrawVTK renderer gcb ${LINK_LIBS} )
endif (WITH_MPI)





###
file (COPY	${CMAKE_CURRENT_SOURCE_DIR}/README.txt
	DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

