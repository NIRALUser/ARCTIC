PROJECT(SegPostProcessCLP)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

GENERATECLP(SegPostProcessCLP.cxx SegPostProcessCLP.xml)
ADD_EXECUTABLE(SegPostProcessCLP SegPostProcessCLP.cxx)
TARGET_LINK_LIBRARIES(SegPostProcessCLP ${ITK_LIBRARIES})

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
