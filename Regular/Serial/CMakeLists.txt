# The sample code to compute entropy for a vector field
cmake_minimum_required (VERSION 2.8)
project (regular_serial_examples) 

find_package(OpenGL)
find_package(GLUT)

# specify the path to the headers
include_directories(
		. 
)

# specify the path to the libraries
link_directories(
)

list(APPEND SRC_FILES
	StreamlineSimple 
	StreamlineSimplePartition 
	)

if( ${WITH_GRAPHICS} )
list(APPEND SRC_FILES
	glStreamlineSimple 
	)
endif()

foreach(SAMPLE_NAME ${SRC_FILES}
	)

	# build the project as a library
	add_executable(${SAMPLE_NAME} ${SAMPLE_NAME}.C)

	if( ${WITH_GRAPHICS} )
		# link the libraries
		if(DEFINED WIN32)
			target_link_libraries(${SAMPLE_NAME} 
				opengl32
				glu32
				glut32
			)
		else()
			target_link_libraries(${SAMPLE_NAME} 
				GL
				GLU
				glut
			)
		endif()
	endif()
	target_link_libraries(${SAMPLE_NAME} 
		OSUFlow
	) 

      # ADD-BY-LEETEN 12/28/2011-BEGIN
      if( ${WITH_PNETCDF} )
	target_link_libraries(${SAMPLE_NAME} 
	  pnetcdf
	  ) 
      endif()

      if( WITH_MPI )
	target_link_libraries(${SAMPLE_NAME} 
	  mpich
	  ) 
      endif( )
      # ADD-BY-LEETEN 12/28/2011-END

endforeach(SAMPLE_NAME)

#
# $Log$
#

