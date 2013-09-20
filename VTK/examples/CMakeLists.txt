# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (VTK_examples) 


set(VTK_DATA_ROOT      "$ENV{HOME}/project/VTKData" CACHE STRING "Root to VTKData for testing the examples")
set(PARAVIEW_DATA_ROOT      "$ENV{HOME}/project/VTKData" CACHE STRING "Root to VTKData for testing the examples")


# specify the path to the headers
include_directories(
	${OSUFLOWVTK_DIR}
)

link_directories (
	${OSUFLOWVTK_DIR}
)

add_definitions (-DVTK_DATA_ROOT=${VTK_DATA_ROOT})

# List all simple targets
LIST(APPEND SIMPLE_TARGETS
	StreamlineSimpleVTK
	testVTS
	vec2vti
	vti2vec
	StreamlineSimpleVTKUnstructured
)
if( ${WITH_GRAPHICS} )
        LIST (APPEND SIMPLE_TARGETS
                glStreamlineSimpleVTK
                gcbStreamlineVTK
                StreamlineSimpleVTKPipeline
                StreamlineRenderVec
                StreamlineRenderPlot3D
		StreamlineRenderUnstructured
        )
endif ()
if ( false ) #${WITH_MPI} )
	LIST (APPEND SIMPLE_TARGETS
        	MpiDrawVTK
	)
	# MPI Script
        file (COPY      ${CMAKE_CURRENT_SOURCE_DIR}/MPI_TEST
                DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
endif ()


# all targets
foreach(TARGET ${SIMPLE_TARGETS} )
	add_executable (${TARGET} ${TARGET}.cpp)
	target_link_libraries (${TARGET} ${LINK_LIBS})
endforeach(TARGET)

########################
# plugin examples
if ( ${WITH_PARAVIEW_PLUGINS} AND ${WITH_MPI})
	include_directories ( ../plugins )
	add_executable (MPIStreamlineRender 
		MPIStreamlineRender.cpp ../plugins/vtkPOSUFlow.cxx)
	target_link_libraries( MPIStreamlineRender ${LINK_LIBS})
endif ()

