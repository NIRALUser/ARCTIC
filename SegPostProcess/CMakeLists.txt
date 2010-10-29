PROJECT(SegPostProcessCLP)

#NeuroLib_USE_ITK()
## NeuroLib_USE_SHAPE()

#FIND_PACKAGE(ITK REQUIRED)
#FIND_PACKAGE(GenerateCLP REQUIRED)

#IF (ITK_FOUND)
#   INCLUDE(${USE_ITK_FILE})
#ELSE(ITK_FOUND)
#   MESSAGE(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
#ENDIF (ITK_FOUND)

#IF (GenerateCLP_FOUND)
#   INCLUDE(${GenerateCLP_USE_FILE})
#ELSE(GenerateCLP_FOUND)
#   MESSAGE(FATAL_ERROR, "GenerateCLP not found. Please set GenerateCLP_DIR.")
#ENDIF(GenerateCLP_FOUND)

SET ( SegPostProcessCLP_SOURCE SegPostProcessCLP.cxx )
GENERATECLP(SegPostProcessCLP_SOURCE SegPostProcessCLP.xml)
ADD_EXECUTABLE ( SegPostProcessCLP ${SegPostProcessCLP_SOURCE})
TARGET_LINK_LIBRARIES ( SegPostProcessCLP ${ITK_LIBRARIES})


slicer3_set_plugins_output_path(SegPostProcessCLP)
slicer3_install_plugins(SegPostProcessCLP)