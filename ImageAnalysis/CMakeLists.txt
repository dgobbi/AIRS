CMAKE_MINIMUM_REQUIRED(VERSION 1.6)

# Set the project/library name here.
PROJECT(vtkImageCompositing)

# Set your list of sources here.  Do not change the name of the
# vtkImageCompsiting_SRCS variable.
SET(vtkImageCompsiting_SRCS
vtkImageCompositeFilter.cxx
)

# List the kits from VTK that are needed by this project.
SET(vtkImageCompsiting_LIBS
  vtkCommon vtkImaging
)

#-----------------------------------------------------------------------------
# Most users should not need to change anything below this line.

# Need to include class headers and the configuration header.
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

# Find and load VTK settings.
INCLUDE(${CMAKE_ROOT}/Modules/FindVTK.cmake)
IF(VTK_FOUND)
  INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
  MESSAGE(FATAL_ERROR "Set VTK_DIR.")
ENDIF(VTK_FOUND)

SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH
     "Single output directory for building all libraries.")

SET (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH
     "Single output directory for building all executables.")

MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

IF(VTK_FOUND)
  # If this is a build tree, provide an option for putting this
  # project's executables and libraries in with VTK's.
  IF (EXISTS ${VTK_DIR}/bin)
    OPTION(USE_VTK_OUTPUT_PATHS
     "Use VTK's output directory for this project's executables and libraries."
     OFF)
    MARK_AS_ADVANCED (USE_VTK_OUTPUT_PATHS)
    IF (USE_VTK_OUTPUT_PATHS)
      SET (LIBRARY_OUTPUT_PATH ${VTK_DIR}/bin)
      SET (EXECUTABLE_OUTPUT_PATH ${VTK_DIR}/bin)
    ENDIF (USE_VTK_OUTPUT_PATHS)
  ENDIF (EXISTS ${VTK_DIR}/bin)

  # Create the library.
  ADD_LIBRARY(${PROJECT_NAME} ${vtkImageCompsiting_SRCS} ${vtkImageCompsitingInstantiator_SRCS})
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${vtkImageCompsiting_LIBS})
  
  # Create Python wrappers if VTK is wrapped in Python.
  IF(VTK_WRAP_PYTHON)
    VTK_WRAP_PYTHON2(${PROJECT_NAME}Python vtkImageCompsitingPYTHON_SRCS ${vtkImageCompsiting_SRCS})
    IF (APPLE)
      ADD_LIBRARY(${PROJECT_NAME}PythonD ${vtkImageCompsitingPYTHON_SRCS})
      ADD_LIBRARY(${PROJECT_NAME}Python MODULE ${PROJECT_NAME}PythonInit.cxx)
      TARGET_LINK_LIBRARIES(${PROJECT_NAME}PythonD ${PROJECT_NAME})
      FOREACH(c ${vtkImageCompsiting_LIBS})
        TARGET_LINK_LIBRARIES(${PROJECT_NAME}PythonD ${c}PythonD)
      ENDFOREACH(c)
      TARGET_LINK_LIBRARIES(${PROJECT_NAME}Python ${PROJECT_NAME}PythonD)
    ELSE (APPLE)
      ADD_LIBRARY(${PROJECT_NAME}Python SHARED ${vtkImageCompsitingPYTHON_SRCS})
      TARGET_LINK_LIBRARIES (${PROJECT_NAME}Python ${PROJECT_NAME})
      FOREACH(c ${vtkImageCompsiting_LIBS})
        TARGET_LINK_LIBRARIES(${PROJECT_NAME}Python ${c}Python)
      ENDFOREACH(c)
    ENDIF(APPLE)

  ENDIF(VTK_WRAP_PYTHON)

  # SUBDIRS(Testing)
ENDIF(VTK_FOUND)
