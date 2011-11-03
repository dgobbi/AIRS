#PROJECT(ITKSegmentation)

SET(KIT ITKSegmentation)
SET(UKIT ITKSEGMENTATION)
SET(KIT_LIBS vtkHybrid)
SET(KIT_EXTRA_LIBS
  ITKIO
  ITKNumerics
  ITKAlgorithms
  ITKCommon
  ITKBasicFilters
  ITKStatistics)

# Set your list of sources here.  Do not change the name of the
# Kit_SRCS variable.
SET(Kit_SRCS
  vtkITKImageAlgorithm.cxx
  vtkITKKnowledgeBasedClassification.cxx
)

SET(Kit_EXTRA_SRCS)
SET(Kit_EXTRA_CMDS)
SET(Kit_TCL_EXTRA_SRCS)
SET(Kit_PYTHON_EXTRA_SRCS)
SET(Kit_JAVA_EXTRA_SRCS)

# -----------------------------------------------------------------
# Include the common setting file
INCLUDE(${AIRS_SOURCE_DIR}/CMake/AIRSCommonBlock.cmake)
