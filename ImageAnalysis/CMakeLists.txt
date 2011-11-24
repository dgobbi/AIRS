#PROJECT(IMAGEANALYSIS)
SET(KIT ImageAnalysis)
SET(UKIT IMAGEANALYSIS)
SET(KIT_LIBS vtkImaging)

# Set your list of sources here.  Do not change the name of the
# Kit_SRCS variable.
SET(Kit_SRCS
vtkImageHistogram.cxx
vtkImageHistogramStatistics.cxx
vtkImageAutoRange.cxx
vtkImageRangeCalculator.cxx
vtkSimpleShapeStencilSource.cxx
)

SET(Kit_EXTRA_SRCS)
SET(Kit_EXTRA_CMDS)
SET(Kit_TCL_EXTRA_SRCS)
SET(Kit_PYTHON_EXTRA_SRCS)
SET(Kit_JAVA_EXTRA_SRCS)

# -----------------------------------------------------------------

# Include the common setting file
INCLUDE(${AIRS_SOURCE_DIR}/CMake/AIRSCommonBlock.cmake)
