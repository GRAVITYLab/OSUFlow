# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (OSUFlow) 

# DEL-BY-LEETEN 08/17/2011	set(MPI FALSE CACHE BOOL "Specify whether MPI is used")

# DEL-BY-LEETEN 01/07/2012-BEGIN
# # MOD-BY-LEETEN 08/17/2011-FROM:
# 	# if( MPI )
# # TO:
# if( WITH_MPI )
# # MOD-BY-LEETEN 08/17/2011-END
# 	# find the mpi package
# 	# DEL-BY-LEETEN 08/17/2011		find_package(MPI REQUIRED)
# 	set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
# 	set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
# 	# DEL-BY-LEETEN 08/17/2011		include_directories(MPI_INCLUDE_PATH)

# 	# use mpicxx to replace 
# 	SET(CMAKE_C_COMPILER mpicc)
# 	SET(CMAKE_CXX_COMPILER mpicxx)
# endif( )
# DEL-BY-LEETEN 01/07/2012-BEGIN

# DEL-BY-LEETEN 08/17/2011-BEGIN
	# # specify the path to the headers
	# include_directories(
		# . 
		# ../include
	# )
# DEL-BY-LEETEN 08/17/2011-END

# grab all the header files
file(GLOB HEADER_FILES "*.h" )

list(APPEND SRC_FILES
	Grid.C       
	polynomials.C  
	TimeVaryingFieldLine.C 
	eigenvals.C  
	Interpolator.C  
	Rake.C	    
	Topology.C 
	eigenvecs.C  
	IsoSurf.C
	Solution.C
	triangulator.C 
	Element.C
	StreakLine.C
	VectorMatrix.C
	Field.C 
	PathLine.C
	Streamline.C
	FieldLine.C
	Plot3DReader.C 
	TimeLine.C
	OSUFlow.C
	FileReader.C
	calc_subvolume.C 
	FlashAMR.C
	ParFlow.C
	flashhdf5_float.C 
	Repartition.C
	ComputeFieldLines.C 
	Blocks.C # ADD-BY-LEETEN 01/06/2012
	Partition.C
	Lattice4D.C
	LatticeAMR.C 
	Plot3DReader.C
	CurvilinearGrid.C
	)

# DEL-BY-LEETEN 01/06/2012-BEGIN
	# # ADD-BY-LEETEN 08/17/2011-BEGIN
	if( WITH_MPI )
	  list(APPEND SRC_FILES
		Blocks.C
	  )
	endif()
	# # ADD-BY-LEETEN 08/17/2011-END
# DEL-BY-LEETEN 01/06/2012-END

# build the project as a library
add_library(${PROJECT_NAME} ${LIB_TYPE} ${SRC_FILES} ${HEADER_FILES})

# ADD-BY-LEETEN 09/01/2011-BEGIN
set_target_properties(${PROJECT_NAME} PROPERTIES 
	DEBUG_OUTPUT_NAME "${PROJECT_NAME}_d"
	RELEASE_OUTPUT_NAME "${PROJECT_NAME}_r"
	)

install(
	FILES ${HEADER_FILES}
	DESTINATION include)
	
install(
	TARGETS ${PROJECT_NAME}
	DESTINATION lib)
# ADD-BY-LEETEN 09/01/2011-END

#
# $Log$
#

