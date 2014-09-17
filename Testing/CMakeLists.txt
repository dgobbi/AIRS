include_directories(${AIRS_INCLUDE_DIRS})
set(AIRS_TESTING_TEMP_DIR ${AIRS_BINARY_DIR}/Testing/Temporary)
set(VTK_TESTING_DIRECTORY "${VTK_DIR}/ExternalData/Testing/")

if(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  set(VTK_LIBS vtkRendering vtkIO)
else(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  set(VTK_LIBS vtkRenderingCore vtkRenderingImage vtkInteractionStyle
      vtkRenderingOpenGL vtkRenderingFreeTypeOpenGL vtkIOImage vtkIOMINC
      vtkImagingStencil vtkImagingMath vtkImagingMorphological
      vtkFiltersSources vtksys)
  # If vtkIOMPIImage is present, it has factories for vtkIOImage
  list(FIND VTK_LIBRARIES vtkIOMPIImage TMP_INDEX)
  if(TMP_INDEX GREATER -1)
    set(VTK_LIBS ${VTK_LIBS} vtkIOMPIImage)
  endif(TMP_INDEX GREATER -1)
endif(${VTK_MAJOR_VERSION} VERSION_LESS 6)

add_executable(TestImageConnectivityFilter
  TestImageConnectivityFilter.cxx)
target_link_libraries(TestImageConnectivityFilter
  vtkImageSegmentation ${VTK_LIBS})
add_test(TestImageConnectivityFilter
  ${CXX_TEST_PATH}/TestImageConnectivityFilter
  -D "${VTK_TESTING_DIRECTORY}")

message(STATUS "TEST ${VTK_TESTING_DIRECTORY}")
