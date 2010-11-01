#
# This module is provided as AIRS_USE_FILE by AIRSConfig.cmake.
# It can be INCLUDEd in a project to load the needed compiler and linker
# settings to use AIRS:
#   FIND_PACKAGE(AIRS REQUIRED)
#   INCLUDE(${AIRS_USE_FILE})
#

IF(NOT AIRS_USE_FILE_INCLUDED)
  SET(AIRS_USE_FILE_INCLUDED 1)

  # Add include directories needed to use AIRS.
  INCLUDE_DIRECTORIES(${AIRS_INCLUDE_DIRS})

  # Add link directories needed to use AIRS.
  LINK_DIRECTORIES(${AIRS_LIBRARY_DIRS})

  #
  # VTK
  #
  IF(NOT VTK_DIR)
    # Use AIRS_VTK_DIR or find a new one
    IF(AIRS_VTK_DIR)
      SET(VTK_DIR ${AIRS_VTK_DIR} CACHE PATH "Path to VTK build dir")
      INCLUDE(${VTK_DIR}/VTKConfig.cmake)
    ELSE(AIRS_VTK_DIR)
      FIND_PACKAGE(VTK REQUIRED)
    ENDIF(AIRS_VTK_DIR)
  ELSE(NOT VTK_DIR)
    INCLUDE(${VTK_DIR}/VTKConfig.cmake)
  ENDIF(NOT VTK_DIR)

  # Include the VTK use file
  INCLUDE(${VTK_USE_FILE})

  #
  # ITK (optional)
  #
  IF(AIRS_USE_ITK)
    IF(NOT ITK_DIR)
      # Use AIRS_ITK_DIR or find a new one
      IF(AIRS_ITK_DIR)
        SET(ITK_DIR ${AIRS_ITK_DIR} CACHE PATH "Path to ITK build dir")
        INCLUDE(${ITK_DIR}/ITKConfig.cmake)
      ELSE(AIRS_ITK_DIR)
        FIND_PACKAGE(ITK REQUIRED)
      ENDIF(AIRS_ITK_DIR)
    ELSE(NOT ITK_DIR)
      INCLUDE(${ITK_DIR}/ITKConfig.cmake)
    ENDIF(NOT ITK_DIR)

    # Include the ITK use file
    INCLUDE(${ITK_USE_FILE})

  ENDIF(AIRS_USE_ITK)

ENDIF(NOT AIRS_USE_FILE_INCLUDED)
