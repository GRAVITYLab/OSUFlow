# The sample code to compute entropy for a vector field
cmake_minimum_required (VERSION 2.8)
project (regular_serial_examples) 

find_package(OpenGL)
find_package(GLUT)

# specify the path to the headers
include_directories(
		. 
		# ADD-BY-LEETEN 01/06/2012-BEGIN
		${RENDERER_DIR}
		${GCB_DIR}
		# ADD-BY-LEETEN 01/06/2012-END
)

# specify the path to the libraries
link_directories(
		# ADD-BY-LEETEN 01/06/2012-BEGIN
		${RENDERER_DIR}
		${GCB_DIR}
		# ADD-BY-LEETEN 01/06/2012-END
)

list(APPEND SRC_FILES
	StreamlineSimple 
	StreamlineSimplePartition 
	PathlineSimple
	# ADD-BY-LEETEN 01/06/2012-BEGIN
	# Incompatible w/ Lattice4D: StreamlineSimplePartition1
	# Incompatible w/ Lattice4D: StreamlineSimplePartition2
	# ADD-BY-LEETEN 01/06/2012-END
	)

if( ${WITH_GRAPHICS} )
list(APPEND SRC_FILES
	glStreamlineSimple 
	# ADD-BY-LEETEN 01/06/2012-BEGIN
	FileDraw
	SerDraw
	glPathlineSimple
	# Incompatible w/ Lattice4D: glPathlinePartition1
	# Incompatible w/ Lattice4D: glPathlinePartition2
	# Incompatible w/ Lattice4D: glPathlinePartitionSimple
	# ADD-BY-LEETEN 01/06/2012-END
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
				# DEL-BY-LEETEN 08/23/2012:	glut32
			)
		else()
			target_link_libraries(${SAMPLE_NAME} 
				GL
				GLU
				glut
			)
		endif()
		# ADD-BY-LEETEN 01/06/2012-BEGIN
		target_link_libraries(${SAMPLE_NAME} 
		  renderer # ADD-BY-LEETEN 01/07/2011
			gcb
			# DEL-BY-LEETEN 01/07/2011: renderer
		)
		# ADD-BY-LEETEN 01/06/2012-END
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
	  # DEL-BY-LEETEN 04/15/2012: mpich
	  diy # ADD-BY-LEETEN 01/07/2011
	  ) 
      endif( )
      # ADD-BY-LEETEN 12/28/2011-END

endforeach(SAMPLE_NAME)

#
# $Log$
#

