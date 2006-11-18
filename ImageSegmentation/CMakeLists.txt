PROJECT(AtamaiSegmentation)

SET(KIT AtamaiSegmentation)
SET(UKIT ATAMAISEGMENTATION)
SET(KIT_LIBS vtkHybrid)

# Set your list of sources here.  Do not change the name of the
# Kit_SRCS variable.
SET(Kit_SRCS
  vtkImageFloodFill.cxx
  #vtkAtamaiPolyDataToImageStencil.cxx Failing on MacOSX-i386, Win32
  vtkAtamaiPolyDataToImageStencil2.cxx
  #vtkImageK9BrainExtractor.cxx - Garbage
)

IF(VTK_MAJOR_VERSION GREATER 4)
  SET(Kit_SRCS ${Kit_SRCS}
    vtkAtamaiMRIBrainExtractor.cxx
    vtkImageMRIBrainExtractor.cxx) # Going to deprecate this!
ELSE(VTK_MAJOR_VERSION GREATER 4)
  SET(Kit_SRCS ${Kit_SRCS}
    vtkImageMRIBrainExtractor.cxx)
ENDIF(VTK_MAJOR_VERSION GREATER 4)

SET(Kit_EXTRA_SRCS)
SET(Kit_EXTRA_CMDS)
SET(Kit_TCL_EXTRA_SRCS)
SET(Kit_PYTHON_EXTRA_SRCS)
SET(Kit_JAVA_EXTRA_SRCS)

# GL: vtkImageMRIBrainExtractor runs ~1.5x slower when compiled with -O3 so,
#     until that gets figured out, force gcc to use -O2.
#SET(CMAKE_CXX_FLAGS_RELEASE "-O2" CACHE STRING "for build with -O2" FORCE)

# -----------------------------------------------------------------

# Include the common setting file
INCLUDE(${AtamaiVTK_SOURCE_DIR}/CMake/KitCommonBlock.cmake)
