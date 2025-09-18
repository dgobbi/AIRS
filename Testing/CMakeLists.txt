include_directories(${AIRS_INCLUDE_DIRS})
set(AIRS_TESTING_TEMP_DIR ${AIRS_BINARY_DIR}/Testing/Temporary)
set(VTK_TESTING_DIRECTORY "${VTK_DIR}/ExternalData/Testing/")

# The library dependencies
set(VTK_TARGET_PREFIX "VTK::")
if(VTK_VERSION VERSION_LESS 8.90)
  set(VTK_TARGET_PREFIX "vtk")
endif()
set(VTK_LIBS)

# required libraries
foreach(_lib RenderingCore RenderingImage InteractionStyle
  RenderingFreeType IOImage IOMINC
  ImagingStencil ImagingMath ImagingMorphological
  FiltersSources)

  list(APPEND VTK_LIBS ${VTK_TARGET_PREFIX}${_lib})
endforeach()

# factory libraries (only needed if present)
foreach(_lib TestingCore IOMPIImage
  RenderingOpenGL RenderingFreeTypeOpenGL
  RenderingOpenGL2 RenderingFreeTypeOpenGL2)

  list(FIND VTK_LIBRARIES ${VTK_TARGET_PREFIX}${_lib} _index)
  if(_index GREATER -1)
    list(APPEND VTK_LIBS ${VTK_TARGET_PREFIX}${_lib})
  endif()
endforeach()

# TODO: This test was moved to VTK proper
add_executable(TestImageConnectivityFilter
  TestImageConnectivityFilter.cxx)
target_link_libraries(TestImageConnectivityFilter
  ${VTK_TARGET_PREFIX}ImageSegmentation ${VTK_LIBS})
get_target_property(_path
  TestImageConnectivityFilter RUNTIME_OUTPUT_DIRECTORY)
add_test(TestImageConnectivityFilter
  "${_path}/TestImageConnectivityFilter"
  -D "${VTK_TESTING_DIRECTORY}")

if(NOT VTK_VERSION VERSION_LESS 8.90)
  vtk_module_autoinit(
    TARGETS TestImageConnectivityFilter
    MODULES VTK::RenderingOpenGL2 VTK::RenderingUI VTK::InteractionStyle)
endif()
