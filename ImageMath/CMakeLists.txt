PROJECT(ImageMath)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

ADD_EXECUTABLE(ImageMath ImageMath.cxx)
TARGET_LINK_LIBRARIES(ImageMath ${ITK_LIBRARIES})

IF (Slicer3_FOUND)
  slicer3_set_plugins_output_path(ImageMath)
  slicer3_install_plugins(ImageMath)
ENDIF (Slicer3_FOUND)

INSTALL(TARGETS
  ImageMath
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static
)
