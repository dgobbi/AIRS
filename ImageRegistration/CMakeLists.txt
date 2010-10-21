#PROJECT(ImageRegistration)

SET(KIT ImageRegistration)
SET(UKIT ImageRegistration)
SET(KIT_LIBS vtkHybrid vtkImaging vtkImageAnalysis)

SET(VTK_WRAP_HINTS ${CMAKE_CURRENT_SOURCE_DIR}/hints)

SET ( Kit_SRCS
vtkCalcCentroid.cxx
vtkCenteredTransform.cxx
vtkImage3DNoiseSource.cxx
vtkImageDataStatistics.cxx
vtkImageMutualInformation.cxx
vtkImageSingleMutualInformation.cxx
vtkCalcCrossCorrelation.cxx
vtkImageAmoebaGrid.cxx
vtkImageGridStatistics.cxx
vtkFunctionMinimizer.cxx
vtkImageCrossCorrelation.cxx
vtkImageMean3D.cxx
)

IF (VTK_USE_ATAMAIIMAGING)
  IF (${VTK_MAJOR_VERSION} GREATER 4)
    SET( Kit_SRCS ${Kit_SRCS}
      vtkImageRegistration.cxx
      )
  ENDIF (${VTK_MAJOR_VERSION} GREATER 4)
ENDIF (VTK_USE_ATAMAIIMAGING )

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
