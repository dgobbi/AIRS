IF (PYTHON_EXECUTABLE)
  ADD_TEST(vtkImageProject-image ${PYTHON_EXECUTABLE}
      ${AtamaiVTK_SOURCE_DIR}/vtkAtamaiImaging/Testing/ImageProject.py
  )
ENDIF (PYTHON_EXECUTABLE)

