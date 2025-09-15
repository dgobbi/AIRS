# Include directories
include_directories(${AIRS_INCLUDE_DIRS})

# Create the main library
set(LIB_NAME vtkImageSegmentation)

# Sources in the current directory (library sources only!)
set(LIB_SRCS
  vtkImageIslandRemoval.cxx
  vtkImageMRIBrainExtractor.cxx
  vtkImageExtractPoints.cxx
  vtkImageRegionIteratorBase.cxx
  vtkImageRegionIterator.cxx
  vtkImagePointsIterator.cxx
)

# Headers without a matching .cxx file are listed here
set(LIB_HDRS
  ${CMAKE_CURRENT_BINARY_DIR}/vtkImageSegmentationModule.h
)

# Internal headers (not installed)
set(LIB_PRIVATE_HDRS
)

# Sources that are abstract
set(LIB_ABSTRACT
  vtkFunctionMinimizer.cxx
  vtkImageSimilarityMetric.cxx
)

# Set source file properties
set_source_files_properties(${LIB_HDRS}
  PROPERTIES WRAP_EXCLUDE ON)
set_source_files_properties(${LIB_ABSTRACT}
  PROPERTIES ABSTRACT ON)
set_source_files_properties(${LIB_PRIVATE_HDRS}
  PROPERTIES SKIP_HEADER_INSTALL ON)

# Create the library
if(DEFINED VTK_MODULE_ENABLE_VTK_ImageSegmentation)

  # When building as a module for VTK 8.90 or later
  set(LIB_SRC_HDRS)
  foreach(_src ${LIB_SRCS})
    get_filename_component(_base ${_src} NAME_WE)
    list(APPEND LIB_SRC_HDRS ${_base}.h)
  endforeach()
  vtk_module_add_module(VTK::ImageSegmentation
    SOURCES ${LIB_SRCS}
    HEADERS ${LIB_SRC_HDRS} ${LIB_HDRS}
    PRIVATE_HEADERS ${LIB_PRIVATE_HDRS})

else()

  # The remainder of the file is for old, non-modular builds
  set(VTK_LIBS
    vtkImagingStencil
    vtkFiltersModeling
    vtkFiltersSources)

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/vtkImageSegmentationModule.h.in
    "${CMAKE_CURRENT_BINARY_DIR}/vtkImageSegmentationModule.h" COPYONLY)

  include("${AIRS_CMAKE_DIR}/AIRSCommonBlock.cmake")

endif()
