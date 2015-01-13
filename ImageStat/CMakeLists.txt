PROJECT(ImageStat)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


ADD_EXECUTABLE( ImageStat ImageStat.cc)
TARGET_LINK_LIBRARIES( ImageStat ${ITK_LIBRARIES} )

IF (Slicer3_FOUND)
  slicer3_set_plugins_output_path(ImageStat)
  slicer3_install_plugins(ImageStat)
ENDIF (Slicer3_FOUND)

INSTALL(TARGETS
  ImageStat
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static
)
