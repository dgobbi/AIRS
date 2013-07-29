PROJECT (Programs)

INCLUDE_DIRECTORIES(${AIRS_INCLUDE_DIRS})

ADD_EXECUTABLE(register register.cxx)
TARGET_LINK_LIBRARIES(register vtkImageRegistration vtkRendering vtkIO ${DICOM_LIBRARIES})
