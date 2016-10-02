#PROJECT(ImageRegistration)

SET(KIT ImageRegistration)
SET(UKIT ImageRegistration)
IF(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  SET(KIT_LIBS vtkHybrid vtkImaging)
ELSE(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  SET(KIT_LIBS vtkImagingCore vtkImagingStatistics vtkImagingStencil)
ENDIF(${VTK_MAJOR_VERSION} VERSION_LESS 6)

SET(VTK_WRAP_HINTS ${CMAKE_CURRENT_SOURCE_DIR}/hints)

SET ( Kit_SRCS
vtkFrameFinder.cxx
vtkFunctionMinimizer.cxx
vtkImageMutualInformation.cxx
vtkImageSquaredDifference.cxx
vtkMorphologicalInterpolator.cxx
vtkLabelInterpolator.cxx
vtkGaussianInterpolator.cxx
vtkImageCorrelationRatio.cxx
vtkImageCrossCorrelation.cxx
vtkImageNeighborhoodCorrelation.cxx
vtkImageSimilarityMetric.cxx
vtkImageRegistrationBase.cxx
vtkImageRegistration.cxx
vtkImageSpread.cxx
vtkITKXFMReader.cxx
vtkITKXFMWriter.cxx
vtkPowellMinimizer.cxx
vtkNelderMeadMinimizer.cxx
)

SET_SOURCE_FILES_PROPERTIES(
vtkFunctionMinimizer
vtkImageRegistrationBase
vtkImageSimilarityMetric
ABSTRACT
)

SET(Kit_EXTRA_SRCS)
SET(Kit_EXTRA_CMDS)
SET(Kit_TCL_EXTRA_SRCS)
SET(Kit_PYTHON_EXTRA_SRCS)
SET(Kit_JAVA_EXTRA_SRCS)

# -----------------------------------------------------------------
# Include the common setting file
INCLUDE(${AIRS_SOURCE_DIR}/CMake/AIRSCommonBlock.cmake)
