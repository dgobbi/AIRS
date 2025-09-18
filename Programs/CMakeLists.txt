include_directories(${AIRS_INCLUDE_DIRS})

# The library dependencies
set(VTK_TARGET_PREFIX "VTK::")
if(VTK_VERSION VERSION_LESS 8.90)
  set(VTK_TARGET_PREFIX "vtk")
endif()
set(VTK_LIBS)

# required libraries
foreach(_lib RenderingCore RenderingImage InteractionStyle
  RenderingFreeType IOImage IOMINC IOGeometry IOLegacy
  ImagingStencil ImagingStatistics CommonSystem CommonCore)

  list(APPEND VTK_LIBS ${VTK_TARGET_PREFIX}${_lib})
endforeach()

# factory libraries (only needed if present)
foreach(_lib IOMPIImage
  RenderingOpenGL RenderingFreeTypeOpenGL
  RenderingOpenGL2 RenderingFreeTypeOpenGL2)

  list(FIND VTK_LIBRARIES ${VTK_TARGET_PREFIX}${_lib} _index)
  if(_index GREATER -1)
    list(APPEND VTK_LIBS ${VTK_TARGET_PREFIX}${_lib})
  endif()
endforeach()

# utility libraries
if(VTK_VERSION VERSION_LESS 8.90)
  list(APPEND VTK_LIBS vtksys)
else()
  list(APPEND VTK_LIBS ${VTK_TARGET_PREFIX}vtksys)
endif()

add_executable(register register.cxx)
target_link_libraries(register ${VTK_TARGET_PREFIX}ImageRegistration
  ${VTK_LIBS} ${DICOM_LIBRARIES} ${GDCM_LIBS})

add_executable(skullstrip skullstrip.cxx)
target_link_libraries(skullstrip ${VTK_TARGET_PREFIX}ImageSegmentation
  ${VTK_LIBS} ${DICOM_LIBRARIES} ${GDCM_LIBS})

install(TARGETS register skullstrip
  EXPORT DICOMTargets
  RUNTIME DESTINATION ${AIRS_RUNTIME_INSTALL_DEST} COMPONENT RuntimeLibraries
  LIBRARY DESTINATION ${AIRS_LIBRARY_INSTALL_DEST} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${AIRS_ARCHIVE_INSTALL_DEST} COMPONENT Development)
