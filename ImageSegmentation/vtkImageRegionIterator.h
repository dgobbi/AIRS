/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageRegionIterator.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageRegionIterator - an image region iterator
// .SECTION Description
// This is an image iterator that can be used to iterate over a
// region of an image.

// .SECTION See also
// vtkImageData vtkImageStencilData vtkImageProgressIterator

#ifndef __vtkImageRegionIterator_h
#define __vtkImageRegionIterator_h

#include "vtkImageRegionIteratorBase.h"

class vtkDataArray;

template<class DType>
class VTK_EXPORT vtkImageRegionIterator : public vtkImageRegionIteratorBase
{
public:
  // Description:
  // Default empty constructor, useful only when creating an array of
  // iterators. Call Initialize() on each iterator before using.
  vtkImageRegionIterator()
    {
    this->PixelIncrement = 0;
    this->BasePointer = 0;
    this->Pointer = 0;
    this->SpanEndPointer = 0;
    }

  // Description:
  // Create an iterator for an extent of an image.  If a stencil is
  // provided, it must have an extent at least as large as the desired
  // extent.
  vtkImageRegionIterator(
    vtkImageData *image, vtkImageStencilData *stencil, const int extent[6],
    vtkAlgorithm *algorithm=0, int threadId=0)
    : vtkImageRegionIteratorBase(image, stencil, extent, algorithm, threadId)
    {
    this->BasePointer = static_cast<DType *>(
      this->GetBasePointer(image, &this->PixelIncrement));
    this->UpdatePointer();
    }

  // Description:
  // Initialize  an iterator for an extent of the image.  If a stencil is
  // provided, it must have an extent at least as large as the desired
  // extent.
  void Initialize(
    vtkImageData *image, vtkImageStencilData *stencil, const int extent[6])
    {
    this->vtkImageRegionIteratorBase::Initialize(image, stencil, extent);
    this->BasePointer = static_cast<DType *>(
      this->GetBasePointer(image, &this->PixelIncrement));
    this->UpdatePointer();
    }

  // Description:
  // Move the iterator to the start of the next span.  A span is a
  // contiguous region over which nothing but the X index changes.
  void NextSpan()
    {
    this->vtkImageRegionIteratorBase::NextSpan();
    this->UpdatePointer();
    }

  // Description:
  // Test if the end of the extent has been reached
  bool IsAtEnd()
    {
    return this->vtkImageRegionIteratorBase::IsAtEnd();
    }

  // Description:
  // Return a pointer to the beginning of the current span.
  DType *BeginSpan()
    {
    return this->Pointer;
    }

  // Description:
  // Return a pointer to the end of the current span.
  DType *EndSpan()
    {
    return this->SpanEndPointer;
    }

protected:

  void UpdatePointer()
    {
    this->Pointer = this->BasePointer + this->Id*this->PixelIncrement;
    this->SpanEndPointer =
      this->BasePointer + this->SpanEnd*this->PixelIncrement;
    }

  // Array information
  int PixelIncrement;

  // Pointers
  DType *BasePointer;       // pointer to the first voxel
  DType *Pointer;           // current iterator position within data
  DType *SpanEndPointer;    // end of current span
};

#endif
// VTK-HeaderTest-Exclude: vtkImageRegionIterator.h
