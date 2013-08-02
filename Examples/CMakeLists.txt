PROJECT (Examples)

# Building outside of AIRS build tree
IF (NOT AIRS_BINARY_DIR)
  MESSAGE(STATUS "Buiding examples outside of AIRS build tree")

  IF (NOT AIRS_DIR)
    FIND_PACKAGE(AIRS REQUIRED)
  ELSE(NOT AIRS_DIR)
    INCLUDE(${AIRS_DIR}/AIRSConfig.cmake)
  ENDIF (NOT AIRS_DIR)

  INCLUDE(${AIRS_USE_FILE})

  # Test output
  SET(AIRS_EXAMPLES_TEMP_DIR ${AIRS_EXAMPLES_BINARY_DIR}/Testing/Temporary)

  # Configure Dart testing support.
  INCLUDE(${CMAKE_ROOT}/Modules/Dart.cmake)

  # Output paths
  SET(EXECUTABLE_OUTPUT_PATH ${AIRS_EXAMPLES_BINARY_DIR}/bin/ CACHE INTERNAL "Single output directory for building all executables.")
  SET(AIRS_EXAMPLES_EXECUTABLE_DIR ${EXECUTABLE_OUTPUT_PATH}/${CMAKE_CFG_INTDIR})
  SET(CXX_TEST_PATH ${AIRS_EXAMPLES_EXECUTABLE_DIR})

# Testing inside the AIRS build tree
ELSE (NOT AIRS_BINARY_DIR)
  MESSAGE(STATUS "Buiding examples inside the AIRS build tree")

  INCLUDE_DIRECTORIES(${AIRS_INCLUDE_DIRS})

  # Test output
  SET(AIRS_EXAMPLES_TEMP_DIR ${AIRS_BINARY_DIR}/Testing/Temporary)
ENDIF (NOT AIRS_BINARY_DIR)

ADD_EXECUTABLE(RigidImageRegistration RigidImageRegistration.cxx)
TARGET_LINK_LIBRARIES(RigidImageRegistration vtkImageRegistration vtkRendering vtkIO ${DICOM_LIBRARIES})

ADD_EXECUTABLE(DifferenceRegistration DifferenceRegistration.cxx)
TARGET_LINK_LIBRARIES(DifferenceRegistration vtkImageRegistration vtkRendering vtkIO ${DICOM_LIBRARIES})

ADD_EXECUTABLE(FrameFinder FrameFinder.cxx)
TARGET_LINK_LIBRARIES(FrameFinder vtkImageRegistration vtkRendering vtkIO)

IF(AIRS_USE_DICOM)
# actually requires NIFTI, which is in DICOM code repository
  ADD_EXECUTABLE(ComputeStrain ComputeStrain.cxx)
  TARGET_LINK_LIBRARIES(ComputeStrain vtkImageRegistration vtkIO ${DICOM_LIBRARIES})
  ADD_EXECUTABLE(WarpImage WarpImage.cxx)
  TARGET_LINK_LIBRARIES(WarpImage vtkImageRegistration vtkIO ${DICOM_LIBRARIES})
ENDIF(AIRS_USE_DICOM)
