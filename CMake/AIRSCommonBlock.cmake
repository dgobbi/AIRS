# Handle common tasks for VTK < 9.0

# Automatically find the header for each cxx file
foreach(arg ${LIB_SRCS})
  get_filename_component(src "${arg}" ABSOLUTE)
  string(REGEX REPLACE "\\.(cxx|c|mm|m)$" ".h" hdr "${src}")
  if("${hdr}" MATCHES "\\.h$" AND EXISTS "${hdr}")
    list(APPEND LIB_HDRS "${hdr}")
  endif()
endforeach()

# Create the hierarchy file
if(BUILD_PYTHON_WRAPPERS)
  set_source_files_properties(${LIB_SPECIAL} PROPERTIES WRAP_SPECIAL ON)
  set(MODULE_HIERARCHY_NAME ${LIB_NAME}Hierarchy)
  # _LINK_DEPENDS is a variable suffix from the VTK module macros
  set(${LIB_NAME}_LINK_DEPENDS ${VTK_LIBS})
  include(vtkWrapHierarchy)
  vtk_wrap_hierarchy(${LIB_NAME} ${CMAKE_CURRENT_BINARY_DIR} "${LIB_SRCS}")
  set(KIT_HIERARCHY_FILE ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_HIERARCHY_NAME}.txt)
  set(LIB_HIERARCHY_STAMP ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_HIERARCHY_NAME}.stamp.txt)
endif()

# Set the library name suffix for VTK libraries
set(LIB_NAME_SUFFIX "-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
if(DEFINED VTK_CUSTOM_LIBRARY_SUFFIX)
  set(LIB_NAME_SUFFIX "${VTK_CUSTOM_LIBRARY_SUFFIX}")
endif()

# Create the library
add_library(${LIB_NAME} ${LIB_SRCS} ${LIB_HIERARCHY_STAMP})
set_target_properties(${LIB_NAME} PROPERTIES
  OUTPUT_NAME ${LIB_NAME}${LIB_NAME_SUFFIX})
set_target_properties(${LIB_NAME} PROPERTIES
  VERSION "${${PROJECT_NAME}_VERSION}"
  SOVERSION "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}")
if(BUILD_PYTHON_WRAPPERS)
  set_target_properties(${LIB_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
target_link_libraries(${LIB_NAME} LINK_PUBLIC ${VTK_LIBS})

# Wrappers
if(BUILD_PYTHON_WRAPPERS)
  set(MODULE_PYTHON_NAME ${LIB_NAME}Python)
  set(${LIB_NAME}_WRAP_DEPENDS ${VTK_LIBS})
  include(vtkWrapPython)
  include_directories(${vtkPython_INCLUDE_DIRS})
  set(VTK_PYTHON_LIBRARIES ${vtkPython_LIBRARIES})
  set(LIB_PYTHON_LIBS vtkWrappingPythonCore)
  vtk_wrap_python3(${MODULE_PYTHON_NAME} LIB_PYTHON_SRCS "${LIB_SRCS}")
  if(AIRS_PYTHON_LIBRARIES)
    # do things the old way, with PythonD libraries
    set(XY) # Get python version, e.g. 38 for python 3.8
    if(vtkPython_LIBRARIES)
      list(GET vtkPython_LIBRARIES 0 TMP_LIB_NAME)
      get_filename_component(TMP_NAME "${TMP_LIB_NAME}" NAME)
      string(REGEX REPLACE "^[^0-9]*([0-9])\\.*([0-9]).*$" "\\1\\2"
        XY "${TMP_NAME}")
      if(NOT XY)
        set(XY)
      endif()
    endif()
    set(LIB_PYTHON_NAME ${LIB_NAME}PythonD)
    set(LIB_PYTHON_OUTPUT_NAME ${LIB_NAME}Python${XY}D)
    set(LIB_PYTHON_OUTPUT_NAME ${LIB_PYTHON_OUTPUT_NAME}${LIB_NAME_SUFFIX})
    foreach(TMP_LIB ${VTK_LIBS})
      set(LIB_PYTHON_LIBS ${LIB_PYTHON_LIBS} ${TMP_LIB}PythonD)
    endforeach()
    add_library(${LIB_PYTHON_NAME} ${LIB_PYTHON_SRCS} ${LIB_PYTHON_EXTRA_SRCS})
    set_target_properties(${LIB_PYTHON_NAME} PROPERTIES
      POSITION_INDEPENDENT_CODE ON)
    set_target_properties(${LIB_PYTHON_NAME} PROPERTIES
      VERSION "${PROJECT_VERSION}"
      SOVERSION "${PROJECT_MAJOR_VERSION}.${PROJECT_MINOR_VERSION}"
      OUTPUT_NAME "${LIB_PYTHON_OUTPUT_NAME}")
    target_link_libraries(${LIB_PYTHON_NAME} LINK_PUBLIC
      ${LIB_NAME} ${LIB_PYTHON_LIBS})
    # On Win32 and Mac, link python library non-private
    if(WIN32 OR APPLE)
      target_link_libraries(${LIB_PYTHON_NAME} LINK_PUBLIC
        ${VTK_PYTHON_LIBRARIES})
    else()
      target_link_libraries(${LIB_PYTHON_NAME} LINK_PRIVATE
        ${VTK_PYTHON_LIBRARIES})
    endif()
    if(KIT_PYTHON_DEPS)
      add_dependencies(${LIB_PYTHON_NAME} ${KIT_PYTHON_DEPS})
    endif()
    add_library(${MODULE_PYTHON_NAME} MODULE ${MODULE_PYTHON_NAME}Init.cxx)
    target_link_libraries(${MODULE_PYTHON_NAME} ${LIB_PYTHON_NAME})
  else()
    # do things the new way, without PythonD libraries
    add_library(${MODULE_PYTHON_NAME} MODULE ${MODULE_PYTHON_NAME}Init.cxx
      ${LIB_PYTHON_SRCS} ${LIB_PYTHON_EXTRA_SRCS})
    target_link_libraries(${MODULE_PYTHON_NAME} LINK_PRIVATE
      ${LIB_NAME} ${LIB_PYTHON_LIBS} ${VTK_PYTHON_LIBRARIES})
    if(KIT_PYTHON_DEPS)
      add_dependencies(${MODULE_PYTHON_NAME} ${KIT_PYTHON_DEPS})
    endif()
  endif()
  set_target_properties(${MODULE_PYTHON_NAME} PROPERTIES PREFIX "")
  if(WIN32 AND NOT CYGWIN)
    set_target_properties(${MODULE_PYTHON_NAME} PROPERTIES SUFFIX ".pyd")
  endif()
  set_target_properties(${MODULE_PYTHON_NAME} PROPERTIES NO_SONAME 1)
endif()

# Set the install rules for the library
install(TARGETS
  ${LIB_NAME} ${LIB_PYTHON_NAME} ${MODULE_PYTHON_NAME}
  EXPORT AIRSTargets
  RUNTIME DESTINATION ${AIRS_RUNTIME_INSTALL_DEST} COMPONENT RuntimeLibraries
  LIBRARY DESTINATION ${AIRS_LIBRARY_INSTALL_DEST} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${AIRS_ARCHIVE_INSTALL_DEST} COMPONENT Development)

install(FILES ${LIB_HDRS}
  DESTINATION ${AIRS_INCLUDE_INSTALL_DEST} COMPONENT Development)
