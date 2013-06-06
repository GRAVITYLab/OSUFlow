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
target_link_libraries ( StreamlineSimpleVTK OSUFlowVTK ${VTK_LIBRARIES} OSUFlow)

add_executable(testVTS testVTS.cpp)
target_link_libraries ( testVTS ${VTK_LIBRARIES} )

add_executable(vec2vts vec2vts.cpp )
target_link_libraries ( vec2vts ${VTK_LIBRARIES} )

add_executable(vec2vti vec2vti.cpp )
target_link_libraries ( vec2vti ${VTK_LIBRARIES} )


if( ${WITH_GRAPHICS} )

	add_executable(glStreamlineSimpleVTK glStreamlineSimpleVTK.cpp)
	target_link_libraries ( glStreamlineSimpleVTK OSUFlowVTK ${VTK_LIBRARIES} OSUFlow GL GLU glut)

	# GCB related application
	include_directories ( ${RENDERER_DIR} ${GCB_DIR} )
	link_directories ( ${RENDERER_DIR} ${GCB_DIR} )

	add_executable(gcbStreamlineVTK gcbStreamlineVTK.cpp)
	target_link_libraries( gcbStreamlineVTK  OSUFlowVTK ${VTK_LIBRARIES} OSUFlow gcb renderer GL GLU glut  )

	add_executable(StreamlineSimpleVTKPipeline StreamlineSimpleVTKPipeline.cpp)
	target_link_libraries( StreamlineSimpleVTKPipeline OSUFlowVTK ${VTK_LIBRARIES} OSUFlow GL GLU glut )

	add_executable(StreamlineRenderVec StreamlineRenderVec.cpp)
	target_link_libraries( StreamlineRenderVec OSUFlowVTK ${VTK_LIBRARIES} OSUFlow GL GLU glut )

	add_executable(StreamlineRenderPlot3D StreamlineRenderPlot3D.cpp)
	target_link_libraries( StreamlineRenderPlot3D OSUFlowVTK ${VTK_LIBRARIES} OSUFlow GL GLU glut )

endif( ${WITH_GRAPHICS} )

### MPI Examples ###
if ( WITH_MPI )
	add_executable(MpiDrawVTK MpiDrawVTK.C)
	target_link_libraries ( MpiDrawVTK OSUFlowVTK ${VTK_LIBRARIES} OSUFlow 
				 diy bil )

	if (WITH_PNETCDF)
		target_link_libraries( MpiDrawVTK pnetcdf )
	endif ()

	if (WITH_GRAPHICS)
		target_link_libraries( MpiDrawVTK renderer gcb )
		
		if(DEFINED WIN32)
		    target_link_libraries(MpiDrawVTK
		        opengl32
		        glu32
		        glut32
		        )
		else()
		      target_link_libraries(MpiDrawVTK
		        GL
		        GLU
		        glut
		        )
		endif()
	endif ()
endif (WITH_MPI)





###
file (COPY	${CMAKE_CURRENT_SOURCE_DIR}/README.txt
	DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
