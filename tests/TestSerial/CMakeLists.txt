# The sample code to compute entropy for a vector field
cmake_minimum_required (VERSION 2.8)
project (TestSerial) 

# specify the path to the headers
include_directories(
		${OSUFLOW_INC}
	)

# build the project as a library
add_executable(TestStreamlineSimple TestStreamlineSimple.C)
add_executable(TestPathlineSimple TestPathlineSimple.C)

# link the libraries
target_link_libraries(TestStreamlineSimple
  OSUFlow
  ) 

target_link_libraries(TestPathlineSimple 
  OSUFlow
  )
