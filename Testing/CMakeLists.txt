include_directories(${AIRS_INCLUDE_DIRS})
set(AIRS_TESTING_TEMP_DIR ${AIRS_BINARY_DIR}/Testing/Temporary)
set(VTK_TESTING_DIRECTORY "${VTK_DIR}/ExternalData/Testing/")

if(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  set(VTK_LIBS vtkRendering vtkIO)
else(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  set(VTK_LIBS vtkRenderingCore vtkRenderingImage vtkInteractionStyle
      vtkRenderingFreeType vtkIOImage vtkIOMINC
      vtkImagingStencil vtkImagingMath vtkImagingMorphological
      vtkFiltersSources vtksys)
  set(VTK_EXTRA_LIBS vtkIOMPIImage
      vtkRenderingOpenGL vtkRenderingFreeTypeOpenGL
      vtkRenderingOpenGL2 vtkRenderingFreeTypeOpenGL2)
  foreach(TMP_LIB ${VTK_EXTRA_LIBS})
    list(FIND VTK_LIBRARIES ${TMP_LIB} TMP_INDEX)
    if(TMP_INDEX GREATER -1)
      set(VTK_LIBS ${VTK_LIBS} ${TMP_LIB})
    endif(TMP_INDEX GREATER -1)
  endforeach(TMP_LIB ${VTK_EXTRA_LIBS})
endif(${VTK_MAJOR_VERSION} VERSION_LESS 6)

add_executable(TestImageConnectivityFilter
  TestImageConnectivityFilter.cxx)
target_link_libraries(TestImageConnectivityFilter
  vtkImageSegmentation ${VTK_LIBS})
add_test(TestImageConnectivityFilter
  ${CXX_TEST_PATH}/TestImageConnectivityFilter
  -D "${VTK_TESTING_DIRECTORY}")
