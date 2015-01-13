PROJECT(itkEMSCLP)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../aniso
  ${CMAKE_CURRENT_SOURCE_DIR}/../basicimg
  ${CMAKE_CURRENT_SOURCE_DIR}/../bias
  ${CMAKE_CURRENT_SOURCE_DIR}/../common
  ${CMAKE_CURRENT_SOURCE_DIR}/../conn
  ${CMAKE_CURRENT_SOURCE_DIR}/../qhull
  ${CMAKE_CURRENT_SOURCE_DIR}/../register
  ${CMAKE_CURRENT_SOURCE_DIR}/../robust
  ${CMAKE_CURRENT_SOURCE_DIR}/../spr
  ${CMAKE_CURRENT_SOURCE_DIR}/../xmlio
)

SET(ITKEMSCLP_SRCS itkEMSCLP.cxx )
GENERATECLP(ITKEMSCLP_SRCS itkEMSCLP.xml )

SET(ITKEMS_SRCS
  ../common/Log.cxx
  ../common/MersenneTwisterRNG.cxx
  ../common/EMSTimer.cxx
  ../common/muFile.cxx
  ../qhull/user.c
  ../qhull/global.c
  ../qhull/stat.c
  ../qhull/io.c
  ../qhull/geom2.c
  ../qhull/poly2.c
  ../qhull/merge.c
  ../qhull/qhull.c
  ../qhull/geom.c
  ../qhull/poly.c
  ../qhull/qset.c
  ../qhull/mem.c
  ../register/AmoebaOptimizer.cxx
  ../register/AtlasRegistrationMethod_float+float.cxx
  ../register/ChainedAffineTransform3D.cxx
  ../register/GradientDescentOptimizer.cxx
  ../register/PairRegistrationMethod_float.cxx
  ../register/PowellOptimizer.cxx
  ../register/SimulatedAnnealingOptimizer.cxx
  ../robust/QHullMSTClusteringProcess.cxx
  ../spr/KMeansEstimator.cxx
  ../xmlio/EMSParametersXMLFile.cxx
  EMSParameters.cxx
  EMSegmentationFilter_float+float.cxx
  filterFloatImages.cxx
  runEMS.cxx
  itkEMSCLPPrimary.cxx
  itkEMSCLPCLP.h
)

##TODO: Need to make a library out of the ITKEMS_SRCS
ADD_LIBRARY(ITKEMSCOMMONLIB ${ITKEMS_SRCS})
TARGET_LINK_LIBRARIES(ITKEMSCOMMONLIB ${ITK_LIBRARIES} )
ADD_EXECUTABLE(itkEMSCLP itkEMSCLP.cxx itkEMSCLPCLP.h)
TARGET_LINK_LIBRARIES(itkEMSCLP ITKEMSCOMMONLIB ${ITK_LIBRARIES} )

IF (Slicer3_FOUND)
  slicer3_set_plugins_output_path(ITKEMSCOMMONLIB)
  slicer3_install_plugins(ITKEMSCOMMONLIB)
  slicer3_set_plugins_output_path(itkEMSCLP)
  slicer3_install_plugins(itkEMSCLP)
ENDIF (Slicer3_FOUND)

