/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageRegionIteratorBase.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageRegionIteratorBase - an image region iterator
// .SECTION Description
// This is an image iterator that can be used to iterate over a
// region of an image.

// .SECTION See also
// vtkImageData vtkImageStencilData vtkImageProgressIterator

#ifndef __vtkImageRegionIteratorBase_h
#define __vtkImageRegionIteratorBase_h

#include "vtkSystemIncludes.h"
class vtkDataArray;
class vtkImageData;
class vtkImageStencilData;
class vtkAlgorithm;

#ifndef __WRAP__
class VTK_EXPORT vtkImageRegionIteratorBase
{
public:
  // Description:
  // Default empty constructor, useful only when creating an array of
  // iterators. Call Initialize() on each iterator before using.
  vtkImageRegionIteratorBase();

  // Description:
  // Create an iterator for an extent of an image.  If a stencil is
  // provided, it must have an extent at least as large as the desired
  // extent.
  vtkImageRegionIteratorBase(
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
  // Get the size of the span.
  int GetSpanSize()
    {
    return (this->SpanEnd - this->Id);
    }

  // Description:
  // Test if the end of the extent has been reached
  bool IsAtEnd()
    {
    return (this->Id == this->End);
    }

  // Description:
  // Get the index at the beginning of the current span.
  void GetIndex(int result[3])
    {
    result[0] = this->Index[0];
    result[1] = this->Index[1];
    result[2] = this->Index[2];
    }

  // Description:
  // Get the index at the beginning of the current span.
  const int *GetIndex()
    {
    return this->Index;
    }

  // Description:
  // The Id at the beginning of the current span.
  vtkIdType GetId()
    {
    return this->Id;
    }

  // Description
  // Get a void pointer and pixel increment for the given Id.
  // The pixel increment is the number of scalar components.
  static void *GetVoidPointer(
    vtkImageData *image, int *pixelIncrement, vtkIdType i);

  // Description
  // Get a void pointer and pixel increment for the given Id.
  // The array should hold point attributes of the image.
  static void *GetVoidPointer(
    vtkDataArray *image, int *pixelIncrement, vtkIdType i);

protected:

  // Description
  // Set all the state variables for the stencil span that includes idX.
  void SetSpanState(int idX);

  // Description
  // Report the progress and do an abort check.  This must be called
  // every time that one row of the image is completed. Only called if
  // Algorithm is not null.
  void ReportProgress();

  vtkIdType  Id;                // the current point Id
  vtkIdType  SpanEnd;           // end of current span
  vtkIdType  RowEnd;            // end of current row
  vtkIdType  SliceEnd;          // end of current slice
  vtkIdType  End;               // end of data

  // Increments
  vtkIdType  RowIncrement;      // to same position in next row
  vtkIdType  SliceIncrement;    // to same position in next slice
  vtkIdType  RowEndIncrement;   // from end of row to start of next row
  vtkIdType  SliceEndIncrement; // from end of slice to start of next slice

  // The extent, adjusted for the stencil
  int        Extent[6];

  // Index-related items
  int        Index[3];
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
#endif

#endif
// VTK-HeaderTest-Exclude: vtkImageRegionIteratorBase.h
