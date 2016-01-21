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

#include "vtkSystemIncludes.h"
class vtkImageData;
class vtkImageStencilData;
class vtkAlgorithm;

template<class DType>
class VTK_EXPORT vtkImageRegionIterator
{
public:
  // Description:
  // Default empty constructor, useful only when creating an array of
  // iterators. Call Initialize() on each iterator before using.
  vtkImageRegionIterator();

  // Description:
  // Create an iterator for an extent of an image.  If a stencil is
  // provided, it must have an extent at least as large as the desired
  // extent.
  vtkImageRegionIterator(
    vtkImageData *image, vtkImageStencilData *stencil, const int extent[6],
    vtkAlgorithm *algorithm=0, int threadId=0);

  // Description:
  // Initialize  an iterator for an extent of the image.  If a stencil is
  // provided, it must have an extent at least as large as the desired
  // extent.
  void Initialize(
    vtkImageData *image, vtkImageStencilData *stencil, const int extent[6]);

  // Description:
  // Check if the iterator is within the stencilled region.  This
  // is updated when NextSpan() is called.
  bool IsInStencil()
    {
    return this->InStencil;
    }

  // Description:
  // Move the iterator to the start of the next span.  A span is a
  // contiguous region over which nothing but the X index changes.
  void NextSpan();

  // Description:
  // Test if the end of the extent has been reached
  bool IsAtEnd()
    {
    return (this->PointId == this->End);
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

  // Description:
  // The X index at the beginning of the current span.
  int GetIndexX()
    {
    return this->IndexX;
    }

  // Description:
  // The Y index at the beginning of the current span.
  int GetIndexY()
    {
    return this->IndexY;
    }

  // Description:
  // The Z index at the beginning of the current span.
  int GetIndexZ()
    {
    return this->IndexZ;
    }

  // Description:
  // The PointId at the beginning of the current span.
  int GetPointId()
    {
    return this->PointId;
    }

protected:

  // Description
  // Set all the state variables for the stencil span that includes idX.
  void SetSpanState(int idX);

  // Description
  // Report the progress and do an abort check.  This must be called
  // every time that one row of the image is completed. Only called if
  // Algorithm is not null.
  void ReportProgress();

  // Pointers
  DType     *BasePointer;       // pointer to the first voxel
  DType     *Pointer;           // current iterator position within data
  DType     *SpanEndPointer;    // end of current span

  vtkIdType  PointId;           // the current point Id
  vtkIdType  SpanEnd;           // end of current span
  vtkIdType  RowEnd;            // end of current row
  vtkIdType  SliceEnd;          // end of current slice
  vtkIdType  End;               // end of data

  // Increments
  vtkIdType  PixelIncrement;    // to next pixel
  vtkIdType  RowIncrement;      // to same position in next row
  vtkIdType  SliceIncrement;    // to same position in next slice
  vtkIdType  RowEndIncrement;   // from end of row to start of next row
  vtkIdType  SliceEndIncrement; // from end of slice to start of next slice

  // The span to cover
  int        MinX;
  int        MaxX;
  int        MinY;
  int        MaxY;
  int        MinZ;
  int        MaxZ;

  // Index-related items
  int        IndexX;
  int        IndexY;
  int        IndexZ;
  int        StartY;

  // Stencil-related items
  bool       HasStencil;
  bool       InStencil;
  int        SpanSliceEndIncrement;
  int        SpanSliceIncrement;
  int        SpanIndex;
  int       *SpanCountPointer;
  int      **SpanListPointer;

  // Progress-related items
  vtkAlgorithm *Algorithm;
  vtkIdType  Count;
  vtkIdType  Target;
};

#ifdef VTK_NO_EXPLICIT_TEMPLATE_INSTANTIATION
#include "vtkImageRegionIterator.txx"
#endif

#endif
// VTK-HeaderTest-Exclude: vtkImageRegionIterator.h
