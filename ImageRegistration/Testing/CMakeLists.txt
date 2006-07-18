IF (PYTHON_EXECUTABLE)
  ADD_TEST(TestImageRegistration ${PYTHON_EXECUTABLE}
    ${AtamaiVTK_SOURCE_DIR}/vtkAtamaiRegistration/Testing/TestImageRegistration.py)
ENDIF (PYTHON_EXECUTABLE)
