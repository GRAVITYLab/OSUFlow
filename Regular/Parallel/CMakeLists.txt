# The sample code to compute entropy for a vector field
cmake_minimum_required (VERSION 2.8)
project (mpitest) 

if (${WITH_GRAPHICS})
 find_package(OpenGL)
 find_package(GLU)
 find_package(GLUT)
endif (${WITH_GRAPHICS})

# DEL-BY-LEETEN 01/07/2012-BEGIN
# if( WITH_MPI )
# 	# find the mpi package
# 	find_package(MPI REQUIRED)
# 	set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
# 	set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
# 	include_directories(MPI_INCLUDE_PATH)

# 	# use mpicxx to replace 
# 	SET(CMAKE_C_COMPILER mpicc)
# 	SET(CMAKE_CXX_COMPILER mpicxx)
# endif( )
# DEL-BY-LEETEN 01/07/2012-END

# specify the path to the headers
include_directories(
		. 
		${RENDERER_DIR}
	)

# build the project as a library
add_executable(${PROJECT_NAME} SimpleMpiDraw.C)

# link the libraries
target_link_libraries(${PROJECT_NAME} 
  OSUFlow
  diy
  bil
  ) 

if( WITH_PNETCDF )
  # specify the path to the libraries
  target_link_libraries(${PROJECT_NAME} 
	pnetcdf
  )
endif()

# ADD-BY-LEETEN 01/07/2012-BEGIN
if( ${WITH_GRAPHICS} )
   target_link_libraries(${PROJECT_NAME} 
      renderer
      gcb
      )

      # ADD-BY-JIMMY START 5/14/2013
      # link the opengl libraries
      target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARY} ${GLU_LIBRARY} ${GLUT_LIBRARY})
      # ADD-BY-JIMMY END 5/14/2013
endif()
# ADD-BY-LEETEN 01/07/2012-END


#
# $Log$
#

