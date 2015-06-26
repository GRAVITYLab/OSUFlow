# The sample code to compute entropy for a vector field
cmake_minimum_required (VERSION 2.8)
project (examples) 

# DEL-BY-LEETEN 08/17/2011-BEGIN
	# # set paths to OSUFlow-related libraries
	# set(OSUFLOW_DIR 	"${CMAKE_SOURCE_DIR}" 		) # CACHE	PATH "The path to OSUFlow.")
	# set(RENDERER_DIR	"${OSUFLOW_DIR}/renderer" 	) # CACHE	PATH "The path to the renderer of OSUFlow.")
	# set(GCB_DIR 		"${RENDERER_DIR}/libgcb" 	) # CACHE	PATH "The path to gcb (a GLUT-related mini window managemenet library).")
	# if(NOT DEFINED WIN32 AND NOT DEFINED MAC_OSX_OMPI AND NOT DEFINED MAC_OSX_MPICH)
	#	add_definitions(
	#		-DLINUX
	#	)
	# endif()
# DEL-BY-LEETEN 08/17/2011-END

# specify the path to the headers
include_directories(
		. 
		# ${GCB_DIR}
		# ${GLUT_INCLUDE_PATH}
	)

# specify the path to the libraries
link_directories(
	# DEL-BY-LEETEN 08/17/2011		${OSUFLOW_DIR}/src
	# ${GCB_DIR}
)

	foreach(SAMPLE_NAME 
			gcbTemplate 
			gcbStreamline
			gcbStreamlineDynamicColor
			gcbStreamlineColoringPoints
			gcbTube
			gcbLineAnimator
			gcbTubeDynamicColor
			gcbPathline
			gcbFileRenderer # ADD-BY-LEETEN 10/02/2012
			gcbStreamlineGaussian
	)

	# build the project as a library
	add_executable(${SAMPLE_NAME} ${SAMPLE_NAME}.cpp)

        target_link_libraries(${SAMPLE_NAME} 
                renderer # ADD-BY-LEETEN 10/22/2012
                OSUFlow
                # DEL-BY-LEETEN 10/22/2012: renderer
                gcb
        )

	# ADD-BY-JIMMY 5/10/2013-BEGIN
        target_link_libraries(${SAMPLE_NAME} ${OPENGL_LIBRARY} ${GLU_LIBRARY} ${GLUT_LIBRARY})

      # ADD-BY-LEETEN 12/20/2011-BEGIN
      if( ${WITH_PNETCDF} )
	target_link_libraries(${SAMPLE_NAME} 
	  pnetcdf
	  ) 
      endif()
      # ADD-BY-LEETEN 12/20/2011-END
	
	# ADD-BY-LEETEN 08/17/2011-BEGIN
	if( WITH_MPI )
		target_link_libraries(${SAMPLE_NAME} 
			mpich
			diy
		) 
	endif( )
	# ADD-BY-LEETEN 08/17/2011-END

endforeach(SAMPLE_NAME)

#
# $Log$
#

