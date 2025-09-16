/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageRegionIterator.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkImageSegmentationModule.h" // For export macro
#include "vtkImageRegionIterator.h"

#ifndef VTK_NO_EXPLICIT_TEMPLATE_INSTANTIATION

template class VTKIMAGESEGMENTATION_EXPORT vtkImageRegionIterator<signed char>;
template class VTKIMAGESEGMENTATION_EXPORT vtkImageRegionIterator<char>;
template class VTKIMAGESEGMENTATION_EXPORT vtkImageRegionIterator<int>;
template class VTKIMAGESEGMENTATION_EXPORT vtkImageRegionIterator<long>;
template class VTKIMAGESEGMENTATION_EXPORT vtkImageRegionIterator<short>;
template class VTKIMAGESEGMENTATION_EXPORT vtkImageRegionIterator<float>;
template class VTKIMAGESEGMENTATION_EXPORT vtkImageRegionIterator<double>;
template class VTKIMAGESEGMENTATION_EXPORT vtkImageRegionIterator<unsigned long>;
template class VTKIMAGESEGMENTATION_EXPORT vtkImageRegionIterator<unsigned short>;
template class VTKIMAGESEGMENTATION_EXPORT vtkImageRegionIterator<unsigned char>;
template class VTKIMAGESEGMENTATION_EXPORT vtkImageRegionIterator<unsigned int>;
#if defined(VTK_TYPE_USE_LONG_LONG)
template class VTKIMAGESEGMENTATION_EXPORT vtkImageRegionIterator<long long>;
template class VTKIMAGESEGMENTATION_EXPORT vtkImageRegionIterator<unsigned long long>;
#endif
#if defined(VTK_TYPE_USE___INT64)
template class VTKIMAGESEGMENTATION_EXPORT vtkImageRegionIterator<__int64>;
template class VTKIMAGESEGMENTATION_EXPORT vtkImageRegionIterator<unsigned __int64>;
#endif

#endif
