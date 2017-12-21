PROJECT (Programs)

INCLUDE_DIRECTORIES(${AIRS_INCLUDE_DIRS})

IF(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  SET(VTK_LIBS vtkRendering vtkIO)
ELSE(${VTK_MAJOR_VERSION} VERSION_LESS 6)
  SET(VTK_LIBS vtkRenderingCore vtkRenderingImage vtkInteractionStyle
      vtkRenderingOpenGL vtkRenderingFreeType vtkIOImage vtkIOMINC
      vtkIOLegacy vtkIOGeometry vtkImagingStencil vtkImagingStatistics
      vtksys)
  # If vtkRenderingFreeTypeOpenGL is present, it must be linked
  list(FIND VTK_LIBRARIES vtkRenderingFreeTypeOpenGL TMP_INDEX)
  if(TMP_INDEX GREATER -1)
    set(VTK_LIBS ${VTK_LIBS} vtkRenderingFreeTypeOpenGL)
  endif()
  # If vtkIOMPIImage is present, it has factories for vtkIOImage
  LIST(FIND VTK_LIBRARIES vtkIOMPIImage TMP_INDEX)
  IF(TMP_INDEX GREATER -1)
    SET(VTK_LIBS ${VTK_LIBS} vtkIOMPIImage)
  ENDIF()
  IF(VTK_MAJOR_VERSION VERSION_GREATER 6 OR VTK_MINOR_VERSION VERSION_GREATER 1)
    ADD_DEFINITIONS(-DVTK_HAS_SLAB_SPACING)
  ENDIF(VTK_MINOR_VERSION VERSION_GREATER 1)
ENDIF(${VTK_MAJOR_VERSION} VERSION_LESS 6)

ADD_EXECUTABLE(register register.cxx)
TARGET_LINK_LIBRARIES(register vtkImageRegistration ${VTK_LIBS} ${DICOM_LIBRARIES})

ADD_EXECUTABLE(skullstrip skullstrip.cxx)
TARGET_LINK_LIBRARIES(skullstrip vtkImageSegmentation ${VTK_LIBS} ${DICOM_LIBRARIES})
