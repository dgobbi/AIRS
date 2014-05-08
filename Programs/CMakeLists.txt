PROJECT (Programs)

INCLUDE_DIRECTORIES(${AIRS_INCLUDE_DIRS})

ADD_EXECUTABLE(register register.cxx)
TARGET_LINK_LIBRARIES(register vtkImageRegistration vtkRendering vtkIO ${DICOM_LIBRARIES})

IF(AIRS_USE_DICOM)
# actually requires NIFTI, which is in DICOM code repository
  ADD_EXECUTABLE(lsmtonifti lsmtonifti.cxx)
  TARGET_LINK_LIBRARIES(lsmtonifti vtkImageRegistration vtkIO ${DICOM_LIBRARIES})
  ADD_EXECUTABLE(niftiview niftiview.cxx)
  TARGET_LINK_LIBRARIES(niftiview vtkImageRegistration vtkIO ${DICOM_LIBRARIES})
  ADD_EXECUTABLE(ComputeStrain ComputeStrain.cxx)
  TARGET_LINK_LIBRARIES(ComputeStrain vtkImageRegistration vtkIO ${DICOM_LIBRARIES})
  ADD_EXECUTABLE(WarpImage WarpImage.cxx)
  TARGET_LINK_LIBRARIES(WarpImage vtkImageRegistration vtkIO ${DICOM_LIBRARIES})
  ADD_EXECUTABLE(StrainGlyph StrainGlyph.cxx)
  TARGET_LINK_LIBRARIES(StrainGlyph vtkImageRegistration vtkIO ${DICOM_LIBRARIES})
ENDIF(AIRS_USE_DICOM)
