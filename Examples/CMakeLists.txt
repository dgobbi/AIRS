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

# The library dependencies
set(VTK_TARGET_PREFIX "VTK::")
if(VTK_VERSION VERSION_LESS 8.90)
  set(VTK_TARGET_PREFIX "vtk")
endif()
set(VTK_LIBS)

# required libraries
foreach(_lib RenderingCore RenderingImage InteractionStyle
  RenderingFreeType IOImage IOMINC CommonSystem CommonCore
  ImagingStencil ImagingMath FiltersSources FiltersGeneral)

  list(APPEND VTK_LIBS ${VTK_TARGET_PREFIX}${_lib})
endforeach()

# factory libraries (only needed if present)
foreach(_lib IOMPIImage
  RenderingOpenGL RenderingFreeTypeOpenGL
  RenderingOpenGL2 RenderingFreeTypeOpenGL2)

  list(FIND VTK_LIBRARIES ${VTK_TARGET_PREFIX}${_lib} _index)
  if(_index GREATER -1)
    list(APPEND VTK_LIBS ${VTK_TARGET_PREFIX}${_lib})
  endif()
endforeach()

# utility libraries
if(VTK_VERSION VERSION_LESS 8.90)
  list(APPEND VTK_LIBS vtksys)
else()
  list(APPEND VTK_LIBS ${VTK_TARGET_PREFIX}vtksys)
endif()

add_executable(RigidImageRegistration RigidImageRegistration.cxx)
target_link_libraries(RigidImageRegistration
  ${VTK_TARGET_PREFIX}ImageRegistration
  ${VTK_LIBS} ${DICOM_LIBRARIES} ${GDCM_LIBS})
add_executable(DifferenceRegistration DifferenceRegistration.cxx)
target_link_libraries(DifferenceRegistration
  ${VTK_TARGET_PREFIX}ImageRegistration
  ${VTK_LIBS} ${DICOM_LIBRARIES} ${GDCM_LIBS})

#add_executable(WarpAnts WarpAnts.cxx)
#target_link_libraries(WarpAnts ${VTK_TARGET_PREFIX}ImageRegistration
#  ${VTK_LIBS} ${DICOM_LIBRARIES} ${GDCM_LIBS})

add_executable(FrameFinder FrameFinder.cxx)
target_link_libraries(FrameFinder ${VTK_TARGET_PREFIX}ImageRegistration
  ${VTK_LIBS} ${DICOM_LIBRARIES} ${GDCM_LIBS})
