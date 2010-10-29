PROJECT(ARCTIC)
# Project: ARCTIC (Automatic Regional Cortical ThICkness)

cmake_minimum_required(VERSION 2.6)

IF(CMAKE_COMPILER_2005)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
ENDIF(CMAKE_COMPILER_2005)

IF (${CMAKE_SOURCE_DIR} STREQUAL ${ARCTIC_SOURCE_DIR})
  SET (LIBRARY_OUTPUT_PATH ${ARCTIC_BINARY_DIR}/lib CACHE PATH "Single output directory for building all libraries.")
  SET (EXECUTABLE_OUTPUT_PATH ${ARCTIC_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")
ENDIF (${CMAKE_SOURCE_DIR} STREQUAL ${ARCTIC_SOURCE_DIR})

#-----------------------------------------------------------------------------
# Configure Dart testing support.
OPTION(BUILD_TESTING "Compile Testing." OFF)
IF(BUILD_TESTING)
  ENABLE_TESTING()
  INCLUDE(CTest)
ENDIF(BUILD_TESTING)

FIND_PACKAGE(Slicer3 REQUIRED NO_DEFAULT_PATH)

IF (Slicer3_FOUND)

  INCLUDE(${Slicer3_USE_FILE})
  slicer3_set_default_install_prefix_for_external_projects()
  FIND_PACKAGE(BatchMake REQUIRED)
  INCLUDE(${BatchMake_USE_FILE})

  SET (ARCTIC_SOURCE ARCTIC.cxx ARCTICParameters.cxx ARCTICComputation.cxx)
  GENERATECLP(ARCTIC.cxx ARCTIC.xml)
  ADD_EXECUTABLE ( ARCTIC ${ARCTIC_SOURCE})
  TARGET_LINK_LIBRARIES ( ARCTIC ${ITK_LIBRARIES} ${BatchMake_LIBRARIES})

  ADD_TEST(ARCTICTest1 ${LIBRARY_OUTPUT_PATH}/Slicer3/Plugins/ARCTIC --help)
  ADD_TEST(ARCTICTest2 ${LIBRARY_OUTPUT_PATH}/Slicer3/Plugins/ARCTIC --version)

  slicer3_set_plugins_output_path(ARCTIC)
  slicer3_install_plugins(ARCTIC)

  SUBDIRS(ImageStat)
  SUBDIRS(ImageMath)
  SUBDIRS(SegPostProcess)
  SUBDIRS(CortThick)
  SUBDIRS(itkEMS)

ELSE (Slicer3_FOUND)

MESSAGE(WARNING "Cannot build with a Slicer3 build tree or a Slicer3 installation. Please set Slicer3_DIR. If using a Slicer3 build, Slicer3_DIR needs to be set to the root of the build tree of Slicer3.  If using an installation of Slicer3, Slicer3_DIR needs to be set to <SlicerInstallation>/lib/Slicer3.")

ENDIF (Slicer3_FOUND)

