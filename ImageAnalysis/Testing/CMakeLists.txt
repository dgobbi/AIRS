IF (PYTHON_EXECUTABLE)
  ADD_TEST(vtkImageProject-image ${PYTHON_EXECUTABLE}
      ${AtamaiVTK_SOURCE_DIR}/vtkAtamaiImaging/Testing/ImageProject.py)

  # Test vtkImageHistogramRange
  ADD_TEST(TestImageRangeCalculator ${PYTHON_EXECUTABLE}
      ${AtamaiVTK_SOURCE_DIR}/vtkAtamaiImaging/Testing/ImageRangeCalculator.py)

ENDIF (PYTHON_EXECUTABLE)

