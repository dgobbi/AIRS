# Generate the AIRSConfig.cmake file in the build tree. This doesnot
# configure one for installation. The file tells external projects how to use
# AIRS.

# Help store a literal dollar in a string.  CMake 2.2 allows escaped
# dollars but we have to support CMake 2.0.
SET(DOLLAR "$")

#-----------------------------------------------------------------------------
# Settings for the build tree.

EXPORT_LIBRARY_DEPENDENCIES(
  ${AIRS_BINARY_DIR}/AIRSLibraryDepends.cmake)

# Set the source dir
SET(AIRS_SOURCE_DIR_CONFIG ${AIRS_SOURCE_DIR})

# The library dependencies file.
SET(AIRS_LIBRARY_DEPENDS_FILE
  ${AIRS_BINARY_DIR}/AIRSLibraryDepends.cmake)

INCLUDE(${CMAKE_ROOT}/Modules/CMakeExportBuildSettings.cmake)

CMAKE_EXPORT_BUILD_SETTINGS(
  ${AIRS_BINARY_DIR}/AIRSBuildSettings.cmake)

# The "use" file.
SET(AIRS_USE_FILE_CONFIG
  ${AIRS_BINARY_DIR}/UseAIRS.cmake)

# The build settings file.
SET(AIRS_BUILD_SETTINGS_FILE_CONFIG
  ${AIRS_BINARY_DIR}/AIRSBuildSettings.cmake)

# The library directories.
SET(AIRS_LIBRARY_DIRS_CONFIG ${AIRS_LIBRARY_DIR})

# The kits.
SET(AIRS_KITS_CONFIG ${AIRS_KITS})

# The libraries.
SET(AIRS_LIBRARIES_CONFIG ${AIRS_LIBRARIES})

# The include directories.
SET(AIRS_INCLUDE_DIRS_CONFIG "")
FOREACH(dir ${AIRS_INCLUDE_DIRS})
  SET(AIRS_INCLUDE_DIRS_CONFIG "${AIRS_INCLUDE_DIRS_CONFIG};${dir}")
ENDFOREACH(dir ${AIRS_INCLUDE_DIRS})

# The VTK options.
SET(AIRS_VTK_DIR_CONFIG ${VTK_DIR})

# The library dependencies file.
#IF(NOT AIRS_NO_LIBRARY_DEPENDS)
#  INCLUDE("@AIRS_LIBRARY_DEPENDS_FILE@")
#ENDIF(NOT AIRS_NO_LIBRARY_DEPENDS)

# Configure AIRSConfig.cmake for the build tree.
CONFIGURE_FILE(
  ${AIRS_SOURCE_DIR}/AIRSConfig.cmake.in
  ${AIRS_BINARY_DIR}/AIRSConfig.cmake @ONLY)

# Configure the UseAIRS file
CONFIGURE_FILE(${AIRS_SOURCE_DIR}/UseAIRS.cmake
               ${AIRS_BINARY_DIR}/UseAIRS.cmake COPYONLY)

