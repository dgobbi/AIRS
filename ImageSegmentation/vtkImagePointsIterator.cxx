/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImagePointsIterator.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkImagePointsIterator.h"
#include "vtkImageData.h"

//----------------------------------------------------------------------------
vtkImagePointsIterator::vtkImagePointsIterator()
{
  this->Origin[0] = 0.0;
  this->Origin[1] = 0.0;
  this->Origin[2] = 0.0;
  this->Spacing[0] = 1.0;
  this->Spacing[1] = 1.0;
  this->Spacing[2] = 1.0;
  this->Position[0] = 0.0;
  this->Position[1] = 0.0;
  this->Position[2] = 0.0;
}

//----------------------------------------------------------------------------
vtkImagePointsIterator::vtkImagePointsIterator(
  vtkImageData *image, const int extent[6], vtkImageStencilData *stencil,
  vtkAlgorithm *algorithm, int threadId)
  : vtkImageRegionIteratorBase(image, extent, stencil, algorithm, threadId)
{
  image->GetOrigin(this->Origin);
  image->GetSpacing(this->Spacing);
  this->UpdatePosition();
}

//----------------------------------------------------------------------------
void vtkImagePointsIterator::Initialize(
  vtkImageData *image, const int extent[6], vtkImageStencilData *stencil,
  vtkAlgorithm *algorithm, int threadId)
{
  this->vtkImageRegionIteratorBase::Initialize(
    image, extent, stencil, algorithm, threadId);
  image->GetOrigin(this->Origin);
  image->GetSpacing(this->Spacing);
  this->UpdatePosition();
}
