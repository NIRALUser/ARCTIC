PROJECT(SegPostProcessCLP)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

FIND_PACKAGE(Slicer3 QUIET NO_DEFAULT_PATH)
IF (Slicer3_FOUND)
  INCLUDE(${Slicer3_USE_FILE})
  slicer3_set_default_install_prefix_for_external_projects()
ELSE (Slicer3_FOUND)
  FIND_PACKAGE(ITK REQUIRED)
  IF (ITK_FOUND)
    INCLUDE(${USE_ITK_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
  ENDIF (ITK_FOUND)

  FIND_PACKAGE(GenerateCLP REQUIRED)
  IF(GenerateCLP_FOUND)
    INCLUDE(${GenerateCLP_USE_FILE})
  ELSE(GenerateCLP_FOUND)
    MESSAGE(FATAL_ERROR, "GenerateCLP not found. Please set GenerateCLP_DIR.")
  ENDIF(GenerateCLP_FOUND)
ENDIF (Slicer3_FOUND)

GENERATECLP(SegPostProcessCLP.cxx SegPostProcessCLP.xml)
ADD_EXECUTABLE(SegPostProcessCLP SegPostProcessCLP.cxx)
TARGET_LINK_LIBRARIES(SegPostProcessCLP ITKIO ITKBasicFilters ITKCommon)

IF (Slicer3_FOUND)
  slicer3_set_plugins_output_path(SegPostProcessCLP)
  slicer3_install_plugins(SegPostProcessCLP)
ENDIF (Slicer3_FOUND)

INSTALL(TARGETS
  SegPostProcessCLP
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static
)
