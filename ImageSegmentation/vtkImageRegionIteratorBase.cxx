/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageRegionIteratorBase.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageRegionIteratorBase.h"
#include "vtkImageData.h"
#include "vtkImageStencilData.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkAlgorithm.h"

//----------------------------------------------------------------------------
class vtkImageStencilIteratorFriendship
{
public:

  static int *GetExtentListLengths(vtkImageStencilData *stencil)
    {
    return stencil->ExtentListLengths;
    }

  static int **GetExtentLists(vtkImageStencilData *stencil)
    {
    return stencil->ExtentLists;
    }
};

//----------------------------------------------------------------------------
vtkImageRegionIteratorBase::vtkImageRegionIteratorBase()
{
  this->PointId = 0;
  this->SpanEnd = 0;
  this->RowEnd = 0;
  this->SliceEnd = 0;
  this->End = 0;

  this->RowEndIncrement = 0;
  this->RowIncrement = 0;
  this->SliceEndIncrement = 0;
  this->SliceIncrement = 0;

  this->MinX = 0;
  this->MaxX = 0;
  this->MinY = 0;
  this->MaxY = 0;
  this->MinZ = 0;
  this->MaxZ = 0;

  this->IndexX = 0;
  this->IndexY = 0;
  this->IndexZ = 0;
  this->StartY = 0;

  this->HasStencil = false;
  this->InStencil = false;
  this->SpanSliceEndIncrement = 0;
  this->SpanSliceIncrement = 0;
  this->SpanIndex = 0;
  this->SpanCountPointer = 0;
  this->SpanListPointer = 0;

  this->Algorithm = 0;
  this->Count = 0;
  this->Target = 0;
}

//----------------------------------------------------------------------------
void vtkImageRegionIteratorBase::Initialize(
  vtkImageData *image, vtkImageStencilData *stencil, const int extent[6])
{
  int dataExtent[6];
  image->GetExtent(dataExtent);

  // Compute the increments for marching through the data.
  this->RowIncrement = dataExtent[1] - dataExtent[0] + 1;
  this->SliceIncrement =
    this->RowIncrement*(dataExtent[3] - dataExtent[2] + 1);

  // Compute the span of the image region to be covered.
  int rowSpan = 0;
  int sliceSpan = 0;
  int volumeSpan = 0;
  this->PointId = 0;

  if (extent[1] >= extent[0] &&
      extent[3] >= extent[2] &&
      extent[5] >= extent[4])
    {
    rowSpan = extent[1] - extent[0] + 1;
    sliceSpan = extent[3] - extent[2] + 1;
    volumeSpan = extent[5] - extent[4] + 1;
    this->PointId = (extent[0] - dataExtent[0]) +
      (extent[2] - dataExtent[2])*this->RowIncrement +
      (extent[4] - dataExtent[4])*this->SliceIncrement;
    }
  // Compute the end increments (continous increments).
  this->RowEndIncrement = this->RowIncrement - rowSpan;
  this->SliceEndIncrement = this->RowEndIncrement +
    this->SliceIncrement - this->RowIncrement*sliceSpan;

  // Get the end pointers for row, slice, and volume.
  this->SpanEnd = this->PointId + rowSpan;
  this->RowEnd = this->PointId + rowSpan;
  this->SliceEnd = this->PointId +
    (this->RowIncrement*sliceSpan - this->RowEndIncrement);
  this->End = this->PointId +
    (this->SliceIncrement*volumeSpan - this->SliceEndIncrement);

  // Save the extent (will be adjusted if there is a stencil).
  this->MinX = extent[0];
  this->MaxX = extent[1];
  this->MinY = extent[2];
  this->MaxY = extent[3];
  this->MinZ = extent[4];
  this->MaxZ = extent[5];

  // For keeping track of the current x,y,z index.
  this->IndexX = this->MinX;
  this->IndexY = this->MinY;
  this->IndexZ = this->MinZ;

  // For resetting the Y index after each slice.
  this->StartY = this->IndexY;

  // Code for when a stencil is provided.
  if (stencil)
    {
    this->HasStencil = true;
    this->InStencil = true;

    this->SpanIndex = 0;
    int stencilExtent[6];
    stencil->GetExtent(stencilExtent);

    // The stencil has a YZ array of span lists, we need increments
    // to get to the next Z position in the YZ array.
    this->SpanSliceIncrement = 0;
    this->SpanSliceEndIncrement = 0;

    if (stencilExtent[3] >= stencilExtent[2] &&
        stencilExtent[5] >= stencilExtent[4])
      {
      this->SpanSliceIncrement = stencilExtent[3] - stencilExtent[2] + 1;
      int botOffset = extent[2] - stencilExtent[2];
      if (botOffset >= 0)
        {
        this->SpanSliceEndIncrement += botOffset;
        }
      int topOffset = stencilExtent[3] - extent[3];
      if (topOffset >= 0)
        {
        this->SpanSliceEndIncrement += topOffset;
        }
      }

    // Find the offset to the start position within the YZ array.
    vtkIdType startOffset = 0;

    int yOffset = extent[2] - stencilExtent[2];
    if (yOffset < 0)
      {
      this->MinY = stencilExtent[2];
      }
    else
      {
      startOffset += yOffset;
      }

    if (stencilExtent[3] <= extent[3])
      {
      this->MaxY = stencilExtent[3];
      }

    int zOffset = extent[4] - stencilExtent[4];
    if (zOffset < 0)
      {
      this->MinZ = stencilExtent[4];
      }
    else
      {
      startOffset += zOffset*this->SpanSliceIncrement;
      }

    if (stencilExtent[5] <= extent[5])
      {
      this->MaxZ = stencilExtent[5];
      }

    if (this->MinY <= this->MaxY &&
        this->MinZ <= this->MaxZ)
      {
      this->SpanCountPointer =
        vtkImageStencilIteratorFriendship::GetExtentListLengths(stencil) +
        startOffset;

      this->SpanListPointer =
        vtkImageStencilIteratorFriendship::GetExtentLists(stencil) +
        startOffset;

      // Holds the current position within the span list for the current row
      this->SetSpanState(this->MinX);
      }
    else
      {
      this->SpanCountPointer = 0;
      this->SpanListPointer = 0;
      this->InStencil = false;
      }
    }
  else
    {
    this->HasStencil = false;
    this->InStencil = true;
    this->SpanSliceEndIncrement = 0;
    this->SpanSliceIncrement = 0;
    this->SpanIndex = 0;
    this->SpanCountPointer = 0;
    this->SpanListPointer = 0;
    }
}

//----------------------------------------------------------------------------
vtkImageRegionIteratorBase::vtkImageRegionIteratorBase(
  vtkImageData *image, vtkImageStencilData *stencil, const int extent[6],
  vtkAlgorithm *algorithm, int threadId)
{
  this->Initialize(image, stencil, extent);

  if (algorithm && threadId == 0)
    {
    this->Algorithm = algorithm;
    vtkIdType maxCount = extent[3] - extent[2] + 1;
    maxCount *= extent[5] - extent[4] + 1;
    this->Target = maxCount/50 + 1;
    this->Count = this->Target*50 - (maxCount/this->Target)*this->Target + 1;
    }
  else
    {
    this->Algorithm = 0;
    this->Target = 0;
    this->Count = 0;
    }
}

//----------------------------------------------------------------------------
void vtkImageRegionIteratorBase::SetSpanState(int idX)
{
  // Find the span that includes idX
  bool inStencil = false;
  int *spans = *this->SpanListPointer;
  int n = *this->SpanCountPointer;
  int i;
  for (i = 0; i < n; i++)
    {
    if (spans[i] > idX)
      {
      break;
      }
    inStencil = !inStencil;
    }

  // Set the primary span state variables
  this->SpanIndex = i;
  this->InStencil = inStencil;

  // Clamp the span end to MaxX+1
  int endIdX = this->MaxX + 1;
  if (i < n && spans[i] <= this->MaxX)
    {
    endIdX = spans[i];
    }

  // Compute the pointers for idX and endIdX
  vtkIdType rowStart =
    this->RowEnd - (this->RowIncrement - this->RowEndIncrement);

  this->PointId = rowStart + (idX - this->MinX);
  this->SpanEnd = rowStart + (endIdX - this->MinX);
}

//----------------------------------------------------------------------------
void vtkImageRegionIteratorBase::NextSpan()
{
  if (this->SpanEnd == this->RowEnd)
    {
    int spanIncr = 1;

    if (this->SpanEnd != this->SliceEnd)
      {
      // Move to the next row
      this->PointId = this->RowEnd + this->RowEndIncrement;
      this->RowEnd += this->RowIncrement;
      this->SpanEnd = this->RowEnd;
      this->IndexY++;
      }
    else if (this->SpanEnd != this->End)
      {
      // Move to the next slice
      this->PointId = this->SliceEnd + this->SliceEndIncrement;
      this->SliceEnd += this->SliceIncrement;
      this->RowEnd = this->PointId +
        (this->RowIncrement - this->RowEndIncrement);
      this->SpanEnd = this->RowEnd;
      this->IndexY = this->StartY;
      this->IndexZ++;
      spanIncr += this->SpanSliceEndIncrement;
      }
    else
      {
      // reached End
      this->PointId = this->End;
      return;
      }

    // Start of next row
    this->IndexX = this->MinX;

    if (this->HasStencil)
      {
      if ((this->IndexY >= this->MinY) &&
          (this->IndexY <= this->MaxY) &&
          (this->IndexZ >= this->MinZ) &&
          (this->IndexZ <= this->MaxZ))
        {
        this->SpanCountPointer += spanIncr;
        this->SpanListPointer += spanIncr;
        this->SetSpanState(this->MinX);
        }
      else
        {
        this->InStencil = false;
        }
      }

    if (this->Algorithm)
      {
      this->ReportProgress();
      }
    }
  else
    {
    // Move to the next span in the current row
    this->PointId = this->SpanEnd;
    int spanCount = *this->SpanCountPointer;
    int endIdX = this->MaxX + 1;
    this->IndexX = endIdX;
    if (this->SpanIndex < spanCount)
      {
      int tmpIdX = (*this->SpanListPointer)[this->SpanIndex];
      if (tmpIdX < endIdX)
        {
        this->IndexX = tmpIdX;
        }
      }

    // Get the index to the start of the span after the next
    this->SpanIndex++;
    if (this->SpanIndex < spanCount)
      {
      int tmpIdX = (*this->SpanListPointer)[this->SpanIndex];
      if (tmpIdX < endIdX)
        {
        endIdX = tmpIdX;
        }
      }

    // Compute the end of the span
    this->SpanEnd = this->RowEnd -
      (this->RowIncrement - this->RowEndIncrement) +
      (endIdX - this->MinX);

    // Flip the state
    this->InStencil = !this->InStencil;
    }
}

//----------------------------------------------------------------------------
void *vtkImageRegionIteratorBase::GetBasePointer(
  vtkImageData *image, int *pixelIncrement)
{
  // Get the data array to use.
  vtkDataArray *array = image->GetPointData()->GetScalars();
  *pixelIncrement = array->GetNumberOfComponents();
  return array->GetVoidPointer(0);
}

//----------------------------------------------------------------------------
void vtkImageRegionIteratorBase::ReportProgress()
{
  if (this->Count % this->Target == 0)
    {
    if (this->Algorithm->GetAbortExecute())
      {
      this->PointId = this->End;
      this->SpanEnd = this->End;
      this->RowEnd = this->End;
      this->SliceEnd = this->End;
      }
    else
      {
      this->Algorithm->UpdateProgress(0.02*(this->Count/this->Target));
      }
    }
  this->Count++;
}
