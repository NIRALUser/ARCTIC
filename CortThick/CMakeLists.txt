PROJECT(CortThickCLP)

#FIND_PACKAGE(ITK REQUIRED)
#IF(ITK_FOUND)
#  INCLUDE(${USE_ITK_FILE})
#ELSE(ITK_FOUND)
#  MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
#ENDIF(ITK_FOUND)

#FIND_PACKAGE(VTK REQUIRED)
#IF (VTK_FOUND)
#  INCLUDE (${VTK_USE_FILE})
#ELSE(VTK_FOUND)
#  MESSAGE(FATAL_ERROR, "VTK not found. Please set VTK_DIR.")
#ENDIF (VTK_FOUND)

#IF (GenerateCLP_FOUND)
#   INCLUDE(${GenerateCLP_USE_FILE})
#ELSE(GenerateCLP_FOUND)
#   MESSAGE(FATAL_ERROR, "GenerateCLP not found. Please set GenerateCLP_DIR.")
#ENDIF(GenerateCLP_FOUND)

INCLUDE_DIRECTORIES(
 ${CMAKE_CURRENT_SOURCE_DIR}/Command
 ${CMAKE_CURRENT_SOURCE_DIR}/DistanceCalculate
 ${CMAKE_CURRENT_SOURCE_DIR}/DistanceCalculate/Array1D
 ${CMAKE_CURRENT_SOURCE_DIR}/DistanceCalculate/Array2D
 ${CMAKE_CURRENT_SOURCE_DIR}/DistanceCalculate/Boundary
 ${CMAKE_CURRENT_SOURCE_DIR}/DistanceCalculate/CalculateAngle
 ${CMAKE_CURRENT_SOURCE_DIR}/DistanceCalculate/HashMap
 ${CMAKE_CURRENT_SOURCE_DIR}/DistanceCalculate/Neighborhood
 ${CMAKE_CURRENT_SOURCE_DIR}/DistanceCalculate/WriteVector
 ${CMAKE_CURRENT_SOURCE_DIR}/Gui
 ${CMAKE_CURRENT_SOURCE_DIR}
 ${CMAKE_CURRENT_BINARY_DIR}
)


SET(CXX_FILES 
	Gui/ImageSyl.cxx
	DistanceCalculate/Boundary/boundary.cxx
	DistanceCalculate/toolsdistance.cxx
	DistanceCalculate/CalculateAngle/CalculateAngle.cxx
	DistanceCalculate/HashMap/HashMapDistance.cxx
	DistanceCalculate/Array1D/Tableau1D.cxx
	DistanceCalculate/Neighborhood/Neighborhood.cxx
	DistanceCalculate/Array2D/Tableau2D.cxx
	DistanceCalculate/WriteVector/Vector3D.cxx
	Command/PipelineCommand.cxx
	Command/CommandLineArgumentParser.cxx	
	CortThickCLP.cxx
	)
	
GENERATECLP(CortThickCLP.cxx CortThickCLP.xml )

ADD_EXECUTABLE(CortThickCLP ${CXX_FILES})
TARGET_LINK_LIBRARIES(CortThickCLP ${ITK_LIBRARIES} ${VTK_LIBRARIES} vtkRendering)

slicer3_set_plugins_output_path(CortThickCLP)
slicer3_install_plugins(CortThickCLP)
