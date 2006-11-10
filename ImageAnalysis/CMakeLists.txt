#PROJECT(ATAMAIIMAGING)
SET(KIT AtamaiImaging)
SET(UKIT ATAMAIIMAGING)
SET(KIT_LIBS vtkImaging)

# Set your list of sources here.  Do not change the name of the
# Kit_SRCS variable.
SET(Kit_SRCS
vtkImageCompositeFilter.cxx
vtkImageProjection.cxx
vtkImageReslice2.cxx
vtkImageRangeCalculator.cxx
)

IF (${VTK_MAJOR_VERSION} GREATER 4)
  SET( Kit_SRCS ${Kit_SRCS}
    #vtkImageComposite.cxx
    )
ENDIF (${VTK_MAJOR_VERSION} GREATER 4)

SET(Kit_EXTRA_SRCS)
SET(Kit_EXTRA_CMDS)
SET(Kit_TCL_EXTRA_SRCS)
SET(Kit_PYTHON_EXTRA_SRCS)
SET(Kit_JAVA_EXTRA_SRCS)

# -----------------------------------------------------------------

# Include the common setting file
INCLUDE(${AtamaiVTK_SOURCE_DIR}/CMake/KitCommonBlock.cmake)
