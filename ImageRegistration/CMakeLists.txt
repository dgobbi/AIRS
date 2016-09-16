#PROJECT(ImageRegistration)

SET(KIT ImageRegistration)
SET(UKIT ImageRegistration)
IF(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  SET(KIT_LIBS vtkHybrid vtkImaging)
ELSE(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  SET(KIT_LIBS vtkImagingCore vtkImagingStatistics)
ENDIF(${VTK_MAJOR_VERSION} VERSION_LESS 6)

SET(VTK_WRAP_HINTS ${CMAKE_CURRENT_SOURCE_DIR}/hints)

SET ( Kit_SRCS
vtkFrameFinder.cxx
vtkImageMutualInformation.cxx
vtkImageSquaredDifference.cxx
vtkMorphologicalInterpolator.cxx
vtkLabelInterpolator.cxx
vtkGaussianInterpolator.cxx
vtkImageCorrelationRatio.cxx
vtkImageCrossCorrelation.cxx
vtkImageNeighborhoodCorrelation.cxx
vtkITKXFMReader.cxx
vtkITKXFMWriter.cxx
vtkPowellMinimizer.cxx
)

IF (${VTK_MAJOR_VERSION} GREATER 4)
  SET( Kit_SRCS ${Kit_SRCS}
    vtkImageRegistration.cxx
    )
ENDIF (${VTK_MAJOR_VERSION} GREATER 4)

##IF (${VTK_MAJOR_VERSION} GREATER 4)
##  SET( Kit_SRCS ${Kit_SRCS}
##    vtkImageMean3D.cxx
##    )
##ENDIF (${VTK_MAJOR_VERSION} GREATER 4)

SET(Kit_EXTRA_SRCS)
SET(Kit_EXTRA_CMDS)
SET(Kit_TCL_EXTRA_SRCS)
SET(Kit_PYTHON_EXTRA_SRCS)
SET(Kit_JAVA_EXTRA_SRCS)

# -----------------------------------------------------------------
# Include the common setting file
INCLUDE(${AIRS_SOURCE_DIR}/CMake/AIRSCommonBlock.cmake)
