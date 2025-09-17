cmake_minimum_required(VERSION 3.12)

project(AIRS)

# Allow custom cmake overrides
set(AIRS_CMAKE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/CMake")
set(CMAKE_MODULE_PATH "${AIRS_CMAKE_DIR}" ${CMAKE_MODULE_PATH})

# Project version
set(AIRS_MAJOR_VERSION 1)
set(AIRS_MINOR_VERSION 0)
set(AIRS_PATCH_VERSION 0)
set(AIRS_SHORT_VERSION "${AIRS_MAJOR_VERSION}.${AIRS_MINOR_VERSION}")
set(AIRS_VERSION "${AIRS_SHORT_VERSION}.${AIRS_PATCH_VERSION}")

# Include directories
set(AIRS_INCLUDE_DIRS "${AIRS_BINARY_DIR}")
if(VTK_VERSION VERSION_LESS 8.90)
  list(APPEND AIRS_INCLUDE_DIRS
    "${AIRS_SOURCE_DIR}/ImageRegistration"
    "${AIRS_BINARY_DIR}/ImageRegistration"
    "${AIRS_SOURCE_DIR}/ImageSegmentation"
    "${AIRS_BINARY_DIR}/ImageSegmentation")
endif()

# Find VTK and get VTK configuration information
find_package(VTK REQUIRED)

# Define VTK_VERSION if not already defined
if(NOT DEFINED VTK_VERSION)
  set(VTK_VERSION ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION})
endif()
if(VTK_VERSION VERSION_LESS 8.90)
  include(${VTK_USE_FILE})
endif()
if(VTK_VERSION VERSION_LESS 7.1)
  message(FATAL_ERROR "VTK 7.1 or later is required.")
endif()

# Look for DICOM add-ons to VTK
option(AIRS_USE_DICOM "Use DICOM in AIRS" OFF)
if(AIRS_USE_DICOM)
  find_package(DICOM QUIET)
  if(DICOM_FOUND)
    include(${DICOM_USE_FILE})
  endif()
  # DICOM configured options with defaults
  # set(DICOM_BUILD_SHARED_LIBS "ON")
  # set(DICOM_USE_ITK_GDCM "OFF")
  # set(DICOM_USE_GDCM "ON")
  # set(DICOM_USE_DCMTK "OFF")
  # set(DICOM_USE_SQLITE "OFF")
  # set(DICOM_USE_VTKZLIB "ON")

  if(DICOM_USE_ITK_GDCM)
    # To avoid conflicts, DICOM might use ITK's GDCM
    find_package(ITK REQUIRED)
    if(ITK_FOUND)
      include(${ITK_USE_FILE})
      set(GDCM_LIBS gdcmMSFF gdcmIOD gdcmDSED gdcmDICT)
      set(GDCM_FOUND 1)
    endif()
  elseif(DICOM_USE_GDCM)
    find_package(GDCM)
    if(NOT GDCM_DIR)
      message(FATAL_ERROR "Please set GDCM_DIR.")
    endif()
    if(GDCM_FOUND)
      include(${GDCM_USE_FILE})
      set(GDCM_LIBS gdcmMSFF gdcmIOD gdcmDSED gdcmDICT)
    endif()
  endif()

  set(DCMTK_LIBS)
  if(DICOM_USE_DCMTK)
    find_package(DCMTK)
    if(NOT DCMTK_FOUND)
      message(FATAL_ERROR "DCMTK not found or incomplete.")
    endif()
    # DCMTK headers might require a specific version of C++
    find_file(_dcmtk_config_header "dcmtk/config/osconfig.h" NO_CACHE
      PATHS ${DCMTK_INCLUDE_DIRS})
    if(_dcmtk_config_header)
      file(READ ${_dcmtk_config_header} _dcmtk_config)
      set(_cxx_config_regex "#define HAVE_CXX([0-9][0-9]) 1")
      string(REGEX MATCHALL "${_cxx_config_regex}" _cxx_config_list "${_dcmtk_config}")
      if(_cxx_config_list)
        list(SORT _cxx_config_list)
        list(GET _cxx_config_list -1 _cxx_config)
        string(REGEX REPLACE "${_cxx_config_regex}" "\\1" _cxx_version "${_cxx_config}")
        if("${_cxx_version}" STRGREATER "${CMAKE_CXX_STANDARD}")
          set(CMAKE_CXX_STANDARD ${_cxx_version})
        endif()
      endif()
    endif()
    if(${DCMTK_charls_LIBRARY})
      set(DCMTK_LIBS ${DCMTK_LIBRARIES} ${DCMTK_charls_LIBRARY})
    else()
      set(DCMTK_LIBS ${DCMTK_LIBRARIES})
    endif()
    if(APPLE)
      list(APPEND DCMTK_LIBS iconv)
    endif()
    include_directories(${DCMTK_INCLUDE_DIRS})
  endif()

  set(SQLITE_LIBS)
  if(DICOM_USE_SQLITE)
    set(SQLITE_LIBS sqlite3)
  endif()

endif()

# Standard option
option(BUILD_SHARED_LIBS "Build shared libraries." OFF)

# Use ctest for dashboard submissions
option(BUILD_TESTING "Build the tests" OFF)
if(BUILD_TESTING)
  include(CTest)
endif()

# Which parts of the package to build
option(BUILD_PROGRAMS "Build some useful AIRS tools" ON)
option(BUILD_EXAMPLES "Build the examples" OFF)

# Wrap Python option
if(VTK_WRAP_PYTHON)
  option(BUILD_PYTHON_WRAPPERS "Build python wrappers" ON)
else()
  unset(BUILD_PYTHON_WRAPPERS CACHE)
endif()

# Configuration header
CONFIGURE_FILE(CMake/AIRSConfig.h.in
  "${AIRS_BINARY_DIR}/AIRSConfig.h" @ONLY)

# Defines CMAKE_INSTALL_LIBDIR etcetera
include(GNUInstallDirs)

# Set up our directory structure for output libraries and binaries
# (Note: these are the build locations, not the install locations)
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${AIRS_BINARY_DIR}/bin")
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  if(UNIX)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${AIRS_BINARY_DIR}/lib")
  else()
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${AIRS_BINARY_DIR}/bin")
  endif()
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${AIRS_BINARY_DIR}/lib")
endif()

# Building for VTK 9.x or later
if(NOT VTK_VERSION VERSION_LESS 8.90)
  # Clear VTK_LIBS, linking will be handled by modular build
  set(VTK_LIBS)

  option(AIRS_VERSIONED_INSTALL "Install with versioned names." ON)
  mark_as_advanced(AIRS_VERSIONED_INSTALL)
  set(vtk_version_suffix "")
  set(airs_library_suffix "")
  if(VTK_VERSIONED_INSTALL OR AIRS_VERSIONED_INSTALL)
    set(vtk_version_suffix "-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
    set(airs_library_suffix "${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
  endif()
  if(DEFINED AIRS_CUSTOM_LIBRARY_SUFFIX)
    set(airs_library_suffix "${AIRS_CUSTOM_LIBRARY_SUFFIX}")
  endif()

  # Optional arguments for vtk_module_build()
  set(module_build_extra_args)
  if(NOT VTK_VERSION VERSION_LESS 9.4.0)
    list(APPEND module_build_extra_args
      USE_FILE_SETS ON)
  endif()
  if(NOT "${airs_library_suffix}" STREQUAL "")
    list(APPEND module_build_extra_args
      LIBRARY_NAME_SUFFIX "${airs_library_suffix}")
  endif()

  # Where to put config files
  set(cmake_destination
    "${CMAKE_INSTALL_LIBDIR}/cmake/airs-${AIRS_SHORT_VERSION}")

  set(airs_module_files
    "${CMAKE_CURRENT_SOURCE_DIR}/ImageRegistration/vtk.module"
    "${CMAKE_CURRENT_SOURCE_DIR}/ImageSegmentation/vtk.module")

  vtk_module_scan(
    MODULE_FILES        ${airs_module_files}
    WANT_BY_DEFAULT     ON
    PROVIDES_MODULES    airs_modules
    ENABLE_TESTS        "${BUILD_TESTING}")

  vtk_module_build(
    MODULES             ${airs_modules}
    INSTALL_EXPORT      AIRS
    ARCHIVE_DESTINATION "lib"
    HEADERS_DESTINATION "include/vtk${vtk_version_suffix}"
    CMAKE_DESTINATION   "${cmake_destination}"
    LICENSE_DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/airs-${AIRS_SHORT_VERSION}"
    HIERARCHY_DESTINATION "${CMAKE_INSTALL_LIBDIR}/vtk${vtk_version_suffix}/hierarchy/AIRS"
    VERSION             "${AIRS_VERSION}"
    SOVERSION           "1"
    ${module_build_extra_args})

  # Disable deprecation warnings for MSVC
  if(MSVC)
    vtk_module_definitions(${airs_modules}
      PRIVATE -D_CRT_SECURE_NO_DEPRECATE)
  endif()

  # Set targets to be exported by airs-config.cmake
  set(AIRS_LINK_TARGETS ${airs_modules})

  if(BUILD_PYTHON_WRAPPERS)
    find_package(Python${VTK_PYTHON_VERSION} QUIET COMPONENTS Development)

    vtk_module_wrap_python(
      MODULES         ${airs_modules}
      TARGET          VTK::airsmodules
      INSTALL_EXPORT  AIRSPython
      PYTHON_PACKAGE  "airs"
      CMAKE_DESTINATION   "${cmake_destination}"
      LIBRARY_DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      SOABI           "${Python${VTK_PYTHON_VERSION}_SOABI}"
      BUILD_STATIC    OFF)

    vtk_module_python_default_destination(python_destination)

    export(
      EXPORT    AIRSPython
      NAMESPACE VTK::
      FILE      "${cmake_destination}/AIRSPython-targets.cmake")
    install(
      EXPORT      AIRSPython
      NAMESPACE   VTK::
      FILE        AIRSPython-targets.cmake
      DESTINATION "${cmake_destination}")
  endif()

  # Generate config files
  include(CMakePackageConfigHelpers)
  configure_package_config_file(CMake/airs-config.cmake.in
    "${cmake_destination}/airs-config.cmake"
    INSTALL_DESTINATION "${cmake_destination}")
  write_basic_package_version_file("${cmake_destination}/airs-config-version.cmake"
    VERSION "${AIRS_MAJOR_VERSION}.${AIRS_MINOR_VERSION}.${AIRS_PATCH_VERSION}"
    COMPATIBILITY AnyNewerVersion)

  # In root of build directory for backwards compatibility
  file(GENERATE
    OUTPUT  "${CMAKE_BINARY_DIR}/airs-config.cmake"
    CONTENT "include(\"${CMAKE_BINARY_DIR}/${cmake_destination}/airs-config.cmake\")\n")
  configure_file(
    "${CMAKE_BINARY_DIR}/${cmake_destination}/airs-config-version.cmake"
    "${CMAKE_BINARY_DIR}/airs-config-version.cmake"
    COPYONLY)

  install(
    FILES       "${CMAKE_BINARY_DIR}/${cmake_destination}/airs-config.cmake"
                "${CMAKE_BINARY_DIR}/${cmake_destination}/airs-config-version.cmake"
    DESTINATION "${cmake_destination}"
    COMPONENT   "development")

endif() #(NOT VTK_VERSION VERSION_LESS 8.90)

if(Module_vtkAIRS)
  set(AIRS_BIN_DIR ${VTK_INSTALL_RUNTIME_DIR})
  set(AIRS_LIB_DIR ${VTK_INSTALL_LIBRARY_DIR})
  set(AIRS_ARC_DIR ${VTK_INSTALL_ARCHIVE_DIR})
  set(AIRS_INC_DIR ${VTK_INSTALL_INCLUDE_DIR})
endif()

# Set up our target directory structure for "make install"
set(AIRS_BIN_DIR "bin") # for executables and ".dll" libraries
if(NOT AIRS_LIB_DIR)
  set(AIRS_LIB_DIR "lib") # for unix/linux/osx shared libraries
endif()
if(NOT AIRS_ARC_DIR)
  set(AIRS_ARC_DIR "${AIRS_LIB_DIR}") # for static or ".lib" libraries
endif()
if(NOT AIRS_INC_DIR)
  set(AIRS_INC_DIR "include") # for header files
endif()
if(NOT AIRS_DATA_DIR)
  set(AIRS_DATA_DIR "share/airs-${AIRS_SHORT_VERSION}")
endif()
if(NOT AIRS_PGK_DIR)
  set(AIRS_PGK_DIR "${AIRS_LIB_DIR}/airs-${AIRS_SHORT_VERSION}/cmake")
endif()

# Set the full paths to the install tree
set(AIRS_RUNTIME_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}/${AIRS_BIN_DIR})
set(AIRS_LIBRARY_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}/${AIRS_LIB_DIR})
set(AIRS_ARCHIVE_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}/${AIRS_ARC_DIR})
set(AIRS_INCLUDE_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}/${AIRS_INC_DIR})
set(AIRS_DATA_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}/${AIRS_DATA_DIR})
set(AIRS_CMAKE_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}/${AIRS_PGK_DIR})

# Libraries
set(AIRS_LIBRARIES
  vtkImageRegistration
  vtkImageSegmentation)

if(VTK_VERSION VERSION_LESS 8.90)
  # External wrapping previous to VTK 9

  set(VTK_WRAP_INCLUDE_DIRS
    ${AIRS_INCLUDE_DIRS}
    ${VTK_WRAP_INCLUDE_DIRS})

  # For the python wrappers
  if(BUILD_PYTHON_WRAPPERS)
    # TODO: better python package destination for VTK 7 and 8
    set(python_destination "lib")
    if(TARGET vtkCommonCorePythonD OR TARGET vtkCommonKitPythonD)
      set(AIRS_PYTHON_LIBRARIES
        vtkImageRegistrationPythonD
        vtkImageSegmentationPythonD)
    endif()
  endif()
endif()

if(BUILD_PYTHON_WRAPPERS)
  string(CONCAT _content
    "from .vtkImageRegistration import *\n\n"
    "from .vtkImageSegmentation import *\n\n"
    "__all__ = ['vtkImageRegistration',\n"
    "           'vtkImageSegmentation']\n"
    "__version__ = \"${AIRS_VERSION}\"\n")
  if(NOT VTK_VERSION VERSION_LESS 8.90)
    string(PREPEND _content "import vtkmodules\n")
  endif()
  file(GENERATE
    OUTPUT  "${CMAKE_BINARY_DIR}/${python_destination}/airs/__init__.py"
    CONTENT "${_content}")
  install(
    FILES       "${CMAKE_BINARY_DIR}/${python_destination}/airs/__init__.py"
    DESTINATION "${python_destination}/airs/")
endif()

# Disable C and C++ deprecation warnings for MSVC 2005 and later
if(MSVC_VERSION EQUAL 1400 OR MSVC_VERSION GREATER 1400 OR MSVC10)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE
                  -D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-D_SCL_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_WARNINGS)
endif()

# Add subdirs if not handled by modular build
if(VTK_VERSION VERSION_LESS 8.90)
  add_subdirectory(ImageRegistration)
  add_subdirectory(ImageSegmentation)
endif()

# Subdirectories
if(BUILD_PROGRAMS)
  add_subdirectory(Programs)
endif()

if(BUILD_EXAMPLES)
  add_subdirectory(Examples)
endif()

if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(Testing)
endif()

# Install files
install(FILES Copyright.txt
  DESTINATION ${AIRS_DATA_INSTALL_DEST})

if(VTK_VERSION VERSION_LESS 8.90)
  #---------------------------------------------------
  # Export information for other projects

  # Add all targets to the build-tree export set
  export(TARGETS ${AIRS_LIBRARIES}
    ${AIRS_PYTHON_LIBRARIES}
    FILE "${AIRS_BINARY_DIR}/AIRSTargets.cmake")

  # Export the package for use from the build-tree
  # (this registers the build-tree with a global CMake-registry)
  export(PACKAGE AIRS)

  # Create the AIRSConfig.cmake and AIRSConfigVersion.cmake files
  # ... for the build tree
  set(INCLUDE_DIRS_CONFIG ${AIRS_INCLUDE_DIRS})
  set(LIBRARY_DIRS_CONFIG "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
  configure_file(CMake/AIRSConfig.cmake.in
    "${AIRS_BINARY_DIR}/AIRSConfig.cmake" @ONLY)
  # ... for the install tree
  file(RELATIVE_PATH REL_INCLUDE_DIR "${AIRS_CMAKE_INSTALL_DEST}"
    "${AIRS_INCLUDE_INSTALL_DEST}")
  set(INCLUDE_DIRS_CONFIG "\${AIRS_PGK_DIR}/${REL_INCLUDE_DIR}")
  file(RELATIVE_PATH REL_LIBRARY_DIR "${AIRS_CMAKE_INSTALL_DEST}"
    "${AIRS_LIBRARY_INSTALL_DEST}")
  set(LIBRARY_DIRS_CONFIG "\${AIRS_PGK_DIR}/${REL_LIBRARY_DIR}")
  configure_file(CMake/AIRSConfig.cmake.in
    "${AIRS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/AIRSConfig.cmake" @ONLY)

  # ... for both
  configure_file(CMake/AIRSConfigVersion.cmake.in
    "${AIRS_BINARY_DIR}/AIRSConfigVersion.cmake" @ONLY)

  # Install the AIRSConfig.cmake and AIRSConfigVersion.cmake
  install(FILES
    "${AIRS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/AIRSConfig.cmake"
    "${AIRS_BINARY_DIR}/AIRSConfigVersion.cmake"
    DESTINATION "${AIRS_CMAKE_INSTALL_DEST}" COMPONENT Development)

  # Install the Use file in the build and install directories
  configure_file(${AIRS_SOURCE_DIR}/CMake/UseAIRS.cmake.in
    ${AIRS_BINARY_DIR}/UseAIRS.cmake COPYONLY IMMEDIATE)
  install(FILES
    "${AIRS_BINARY_DIR}/UseAIRS.cmake"
    DESTINATION "${AIRS_CMAKE_INSTALL_DEST}" COMPONENT Development)

  # Install the export set for use with the install-tree
  install(EXPORT AIRSTargets
    DESTINATION "${AIRS_CMAKE_INSTALL_DEST}" COMPONENT Development)
endif()
