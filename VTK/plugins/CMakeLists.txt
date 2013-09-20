cmake_minimum_required(VERSION 2.8)

# The followings are pre-defined in the upper folder
#FIND_PACKAGE(ParaView REQUIRED)
#INCLUDE(${PARAVIEW_USE_FILE})
#include ( ParaViewBranding )
#include ( ParaViewPlugins )

SET(PLUGIN_NAME OSUFlowVTKPlugins)
SET(PLUGIN_VERSION "0.1")

SET(SM_XML GenStreamlines.xml )

SET(SM_SRC
  ../vtkOSUFlow.cxx
  vtkPOSUFlow.cxx
  )

ADD_PARAVIEW_PLUGIN(${PLUGIN_NAME} ${PLUGIN_VERSION}
  REQUIRED_ON_SERVER
  SERVER_MANAGER_XML ${SM_XML}
  SERVER_MANAGER_SOURCES ${SM_SRC}
  )

TARGET_LINK_LIBRARIES(${PLUGIN_NAME} ${LINK_LIBS})

