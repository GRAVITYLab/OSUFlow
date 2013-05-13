# Install script for directory: /home/jchen/project/osuflow/trunk/src

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/jchen/project/osuflow/trunk/src/header.h"
    "/home/jchen/project/osuflow/trunk/src/Repartition.h"
    "/home/jchen/project/osuflow/trunk/src/Field.h"
    "/home/jchen/project/osuflow/trunk/src/polynomials.h"
    "/home/jchen/project/osuflow/trunk/src/FieldLine.h"
    "/home/jchen/project/osuflow/trunk/src/calc_subvolume.h"
    "/home/jchen/project/osuflow/trunk/src/flashhdf5_float.h"
    "/home/jchen/project/osuflow/trunk/src/Partition.h"
    "/home/jchen/project/osuflow/trunk/src/IsoSurf.h"
    "/home/jchen/project/osuflow/trunk/src/Element.h"
    "/home/jchen/project/osuflow/trunk/src/Rake.h"
    "/home/jchen/project/osuflow/trunk/src/Interpolator.h"
    "/home/jchen/project/osuflow/trunk/src/Topology.h"
    "/home/jchen/project/osuflow/trunk/src/CurvilinearGrid.h"
    "/home/jchen/project/osuflow/trunk/src/CandidateCP.h"
    "/home/jchen/project/osuflow/trunk/src/Solution.h"
    "/home/jchen/project/osuflow/trunk/src/triangulator.h"
    "/home/jchen/project/osuflow/trunk/src/Grid.h"
    "/home/jchen/project/osuflow/trunk/src/mcube.h"
    "/home/jchen/project/osuflow/trunk/src/Lattice.h"
    "/home/jchen/project/osuflow/trunk/src/Plot3DReader.h"
    "/home/jchen/project/osuflow/trunk/src/VectorMatrix.h"
    "/home/jchen/project/osuflow/trunk/src/FlashAMR.h"
    "/home/jchen/project/osuflow/trunk/src/LatticeAMR.h"
    "/home/jchen/project/osuflow/trunk/src/ParFlow.h"
    "/home/jchen/project/osuflow/trunk/src/OSUFlow.h"
    "/home/jchen/project/osuflow/trunk/src/zero_test.h"
    "/home/jchen/project/osuflow/trunk/src/Lattice4D.h"
    "/home/jchen/project/osuflow/trunk/src/Newton.h"
    "/home/jchen/project/osuflow/trunk/src/FileReader.h"
    "/home/jchen/project/osuflow/trunk/src/definition.h"
    "/home/jchen/project/osuflow/trunk/src/ComputeFieldLines.h"
    "/home/jchen/project/osuflow/trunk/src/Blocks.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/jchen/project/osuflow/trunk/src/libOSUFlow.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

