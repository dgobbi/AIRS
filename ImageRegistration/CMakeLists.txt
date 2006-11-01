PROJECT(ATAMAIREGISTRATION)

SET(KIT AtamaiRegistration)
SET(UKIT ATAMAIREGISTRATION)
SET(KIT_LIBS vtkHybrid)

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
vtkImageRegistration.cxx
)

SET(Kit_EXTRA_SRCS)
SET(Kit_EXTRA_CMDS)
SET(Kit_TCL_EXTRA_SRCS)
SET(Kit_PYTHON_EXTRA_SRCS)
SET(Kit_JAVA_EXTRA_SRCS)

# -----------------------------------------------------------------

# Include the common setting file
INCLUDE(${AtamaiVTK_SOURCE_DIR}/CMake/KitCommonBlock.cmake)