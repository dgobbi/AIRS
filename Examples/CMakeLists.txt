project(AIRSExamples)

# Building outside of AIRS build tree
if(NOT AIRS_BINARY_DIR)
  message(STATUS "Buiding examples outside of AIRS build tree")

  IF(NOT AIRS_DIR)
    find_package(AIRS REQUIRED)
  else(NOT AIRS_DIR)
    include(${AIRS_DIR}/AIRSConfig.cmake)
  endif()

  include(${AIRS_USE_FILE})

  # Test output
  set(AIRS_EXAMPLES_TEMP_DIR ${AIRS_EXAMPLES_BINARY_DIR}/Testing/Temporary)

  # Output paths
  set(EXECUTABLE_OUTPUT_PATH ${AIRS_EXAMPLES_BINARY_DIR}/bin/ CACHE INTERNAL "Single output directory for building all executables.")
  set(AIRS_EXAMPLES_EXECUTABLE_DIR ${EXECUTABLE_OUTPUT_PATH}/${CMAKE_CFG_INTDIR})
  set(CXX_TEST_PATH ${AIRS_EXAMPLES_EXECUTABLE_DIR})

# Testing inside the AIRS build tree
else(NOT AIRS_BINARY_DIR)
  message(STATUS "Buiding examples inside the AIRS build tree")

  include_directories(${AIRS_INCLUDE_DIRS})

  # Test output
  set(AIRS_EXAMPLES_TEMP_DIR ${AIRS_BINARY_DIR}/Testing/Temporary)
endif()

if(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  set(VTK_LIBS vtkRendering vtkIO)
else(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  set(VTK_LIBS vtkRenderingCore vtkRenderingImage vtkInteractionStyle
      vtkRenderingFreeType vtkIOImage vtkIOMINC vtkCommonSystem
      vtkImagingStencil vtkImagingMath vtkFiltersSources vtkFiltersGeneral
      vtksys)
  set(VTK_EXTRA_LIBS vtkIOMPIImage
      vtkRenderingOpenGL vtkRenderingFreeTypeOpenGL
      vtkRenderingOpenGL2 vtkRenderingFreeTypeOpenGL2)
  foreach(TMP_LIB ${VTK_EXTRA_LIBS})
    list(FIND VTK_LIBRARIES ${TMP_LIB} TMP_INDEX)
    if(TMP_INDEX GREATER -1)
      set(VTK_LIBS ${VTK_LIBS} ${TMP_LIB})
    endif()
  endforeach()
endif()

add_executable(RigidImageRegistration RigidImageRegistration.cxx)
target_link_libraries(RigidImageRegistration
  vtkImageRegistration ${VTK_LIBS} ${DICOM_LIBRARIES})
add_executable(DifferenceRegistration DifferenceRegistration.cxx)
target_link_libraries(DifferenceRegistration
  vtkImageRegistration ${VTK_LIBS} ${DICOM_LIBRARIES})

#add_executable(WarpAnts WarpAnts.cxx)
#target_link_libraries(WarpAnts
#  vtkImageRegistration ${VTK_LIBS} ${DICOM_LIBRARIES})

add_executable(FrameFinder FrameFinder.cxx)
target_link_libraries(FrameFinder
  vtkImageRegistration ${VTK_LIBS} ${DICOM_LIBRARIES})
