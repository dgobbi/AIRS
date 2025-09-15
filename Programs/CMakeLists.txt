include_directories(${AIRS_INCLUDE_DIRS})

set(VTK_LIBS vtkRenderingCore vtkRenderingImage vtkInteractionStyle
    vtkRenderingFreeType vtkIOImage vtkIOMINC vtkIOGeometry vtkIOLegacy
    vtkImagingStencil vtkImagingStatistics vtkCommonSystem vtksys)
set(VTK_EXTRA_LIBS vtkIOMPIImage
    vtkRenderingOpenGL vtkRenderingFreeTypeOpenGL
    vtkRenderingOpenGL2 vtkRenderingFreeTypeOpenGL2)
foreach(TMP_LIB ${VTK_EXTRA_LIBS})
  list(FIND VTK_LIBRARIES ${TMP_LIB} TMP_INDEX)
  if(TMP_INDEX GREATER -1)
    set(VTK_LIBS ${VTK_LIBS} ${TMP_LIB})
  endif()
endforeach()

# TODO: remove corresponding conditionals in register.cxx
add_definitions(-DVTK_HAS_SLAB_SPACING)

add_executable(register register.cxx)
target_link_libraries(register
  vtkImageRegistration ${VTK_LIBS} ${DICOM_LIBRARIES} ${GDCM_LIBS})

add_executable(skullstrip skullstrip.cxx)
target_link_libraries(skullstrip
  vtkImageSegmentation ${VTK_LIBS} ${DICOM_LIBRARIES} ${GDCM_LIBS})

install(TARGETS register skullstrip
  EXPORT DICOMTargets
  RUNTIME DESTINATION ${AIRS_RUNTIME_INSTALL_DEST} COMPONENT RuntimeLibraries
  LIBRARY DESTINATION ${AIRS_LIBRARY_INSTALL_DEST} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${AIRS_ARCHIVE_INSTALL_DEST} COMPONENT Development)
