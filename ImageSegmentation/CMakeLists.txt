#PROJECT(ImageSegmentation)

SET(KIT ImageSegmentation)
SET(UKIT IMAGESEGMENTATION)
IF(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  SET(KIT_LIBS vtkHybrid)
ELSE(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  SET(KIT_LIBS vtkImagingStencil vtkFiltersModeling vtkFiltersSources)
ENDIF(${VTK_MAJOR_VERSION} VERSION_LESS 6)

# Set your list of sources here.  Do not change the name of the
# Kit_SRCS variable.
SET(Kit_SRCS
  vtkImageIslandRemoval.cxx
  vtkImageMRIBrainExtractor.cxx
  vtkImageExtractPoints.cxx
  vtkImageRegionIteratorBase.cxx
  vtkImageRegionIterator.cxx
  vtkImagePointsIterator.cxx
)

SET_SOURCE_FILES_PROPERTIES(
  vtkImageRegionIteratorBase
  vtkImageRegionIterator
  vtkImagePointsIterator
  WRAP_EXCLUDE
)

SET_SOURCE_FILES_PROPERTIES(
  vtkImageRegionIteratorBase
  vtkImageRegionIterator
  vtkImagePointsIterator
  WRAP_PYTHON_EXCLUDE
)

SET(Kit_EXTRA_SRCS)
SET(Kit_EXTRA_CMDS)
SET(Kit_TCL_EXTRA_SRCS)
SET(Kit_PYTHON_EXTRA_SRCS)
SET(Kit_JAVA_EXTRA_SRCS)

# -----------------------------------------------------------------
# Include the common setting file
INCLUDE(${AIRS_SOURCE_DIR}/CMake/AIRSCommonBlock.cmake)
