/*=========================================================================

  Module: vtkImageNeighborhoodCorrelation.cxx

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageNeighborhoodCorrelation.h"

#include "vtkImageSimilarityMetricInternals.h"

#include <vtkObjectFactory.h>
#include <vtkImageData.h>
#include <vtkImageStencilData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTemplateAliasMacro.h>
#include <vtkVersion.h>

// turn off 64-bit ints when templating over all types
# undef VTK_USE_INT64
# define VTK_USE_INT64 0
# undef VTK_USE_UINT64
# define VTK_USE_UINT64 0

#include <math.h>

vtkStandardNewMacro(vtkImageNeighborhoodCorrelation);

//----------------------------------------------------------------------------
// Data needed for each thread.
class vtkImageNeighborhoodCorrelationThreadData
{
public:
  vtkImageNeighborhoodCorrelationThreadData() : Result(0.0) {}

  double Result;
};

class vtkImageNeighborhoodCorrelationTLS
  : public vtkImageSimilarityMetricTLS<vtkImageNeighborhoodCorrelationThreadData>
{
};

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageNeighborhoodCorrelation::vtkImageNeighborhoodCorrelation()
{
  this->NeighborhoodRadius[0] = 7;
  this->NeighborhoodRadius[1] = 7;
  this->NeighborhoodRadius[2] = 7;
}

//----------------------------------------------------------------------------
vtkImageNeighborhoodCorrelation::~vtkImageNeighborhoodCorrelation()
{
}

//----------------------------------------------------------------------------
void vtkImageNeighborhoodCorrelation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "NeighborhoodRadius: " << this->NeighborhoodRadius[0] << " "
     << this->NeighborhoodRadius[1] << " " << this->NeighborhoodRadius[2] << "\n";
}

// begin anonymous namespace
namespace {

//----------------------------------------------------------------------------
// Compute partial sums of x, y, x^2, y*2, x*y for a row of the image,
// given a neighborhood size to use.  Use a sliding-window filter, which
// reduces the computation cost from O(N*M) to O(N) where N is the image
// size and M is the neighborhood size.
template<class T, class U>
void vtkImageNeighborhoodCorrelationX(
  const T *inPtr1, const T *inPtr2, vtkIdType inIncX1, vtkIdType inIncX2,
  int n, int radius, U *workPtr)
{
  if (radius == 0 || n <= 3*radius + 2)
    {
    // use a simple loop with O(n*radius) efficiency
    for (int k = 0; k < n; k++)
      {
      workPtr[0] = 0;
      workPtr[1] = 0;
      workPtr[2] = 0;
      workPtr[3] = 0;
      workPtr[4] = 0;
      workPtr[5] = 0;

      int i = k - radius;
      i = ((i > 0) ? i : 0);
      int j = k + radius + 1;
      j = ((j <= n) ? j : n);
      const T *tmpPtr1 = inPtr1 + i*inIncX1;
      const T *tmpPtr2 = inPtr2 + i*inIncX2;

      for (; i < j; i++)
        {
        U x = *tmpPtr1;
        U y = *tmpPtr2;
        U xx = x*x;
        U yy = y*y;
        U xy = x*y;
        workPtr[0] += x;
        workPtr[1] += y;
        workPtr[2] += xx;
        workPtr[3] += yy;
        workPtr[4] += xy;
        workPtr[5] += 1;
        tmpPtr1 += inIncX1;
        tmpPtr2 += inIncX2;
        }

      workPtr += 6;
      }

    return;
    }

  // use a sliding window to achieve O(n) efficiency rather than O(n*radius)
  workPtr[0] = 0;
  workPtr[1] = 0;
  workPtr[2] = 0;
  workPtr[3] = 0;
  workPtr[4] = 0;
  workPtr[5] = 0;

  // prime
  U *headPtr = workPtr + 6*(radius + 1);
  int k = radius + 1;
  do
    {
    U x = *inPtr1;
    U y = *inPtr2;
    U xx = x*x;
    U yy = y*y;
    U xy = x*y;
    headPtr[0] = x;
    workPtr[0] += x;
    headPtr[1] = y;
    workPtr[1] += y;
    headPtr[2] = xx;
    workPtr[2] += xx;
    headPtr[3] = yy;
    workPtr[3] += yy;
    headPtr[4] = xy;
    workPtr[4] += xy;
    headPtr[5] = 1;
    workPtr[5] += 1;
    headPtr += 6;
    inPtr1 += inIncX1;
    inPtr2 += inIncX2;
    }
  while (--k);

  workPtr += 6;

  // lead
  k = radius;
  do
    {
    U x = *inPtr1;
    U y = *inPtr2;
    U xx = x*x;
    U yy = y*y;
    U xy = x*y;
    headPtr[0] = x;
    workPtr[0] = workPtr[-6] + x;
    headPtr[1] = y;
    workPtr[1] = workPtr[-5] + y;
    headPtr[2] = xx;
    workPtr[2] = workPtr[-4] + xx;
    headPtr[3] = yy;
    workPtr[3] = workPtr[-3] + yy;
    headPtr[4] = xy;
    workPtr[4] = workPtr[-2] + xy;
    headPtr[5] = 1;
    workPtr[5] = workPtr[-1] + 1;
    workPtr += 6;
    headPtr += 6;
    inPtr1 += inIncX1;
    inPtr2 += inIncX2;
    }
  while (--k);

  // slide
  k = n - 3*radius - 2;
  do
    {
    U x = *inPtr1;
    U y = *inPtr2;
    U xx = x*x;
    U yy = y*y;
    U xy = x*y;
    headPtr[0] = x;
    workPtr[0] = workPtr[-6] + x - workPtr[0];
    headPtr[1] = y;
    workPtr[1] = workPtr[-5] + y - workPtr[1];
    headPtr[2] = xx;
    workPtr[2] = workPtr[-4] + xx - workPtr[2];
    headPtr[3] = yy;
    workPtr[3] = workPtr[-3] + yy - workPtr[3];
    headPtr[4] = xy;
    workPtr[4] = workPtr[-2] + xy - workPtr[4];
    headPtr[5] = 1;
    workPtr[5] = workPtr[-1] + 1 - workPtr[5];
    workPtr += 6;
    headPtr += 6;
    inPtr1 += inIncX1;
    inPtr2 += inIncX2;
    }
  while (--k);

  // tail
  k = radius + 1;
  do
    {
    U x = *inPtr1;
    U y = *inPtr2;
    U xx = x*x;
    U yy = y*y;
    U xy = x*y;
    workPtr[0] = workPtr[-6] + x - workPtr[0];
    workPtr[1] = workPtr[-5] + y - workPtr[1];
    workPtr[2] = workPtr[-4] + xx - workPtr[2];
    workPtr[3] = workPtr[-3] + yy - workPtr[3];
    workPtr[4] = workPtr[-2] + xy - workPtr[4];
    workPtr[5] = workPtr[-1] + 1 - workPtr[5];
    workPtr += 6;
    inPtr1 += inIncX1;
    inPtr2 += inIncX2;
    }
  while (--k);

  // finish
  k = radius;
  do
    {
    workPtr[0] = workPtr[-6] - workPtr[0];
    workPtr[1] = workPtr[-5] - workPtr[1];
    workPtr[2] = workPtr[-4] - workPtr[2];
    workPtr[3] = workPtr[-3] - workPtr[3];
    workPtr[4] = workPtr[-2] - workPtr[4];
    workPtr[5] = workPtr[-1] - workPtr[5];
    workPtr += 6;
    }
  while (--k);
}

//----------------------------------------------------------------------------
// Use a stencil to decide which portions of an image row to compute
// the partial sums of x, y, x^2, y^2, and x*y.  The sums are set to
// zero for any voxels outside of the stencil.
template<class T, class U>
void vtkImageNeighborhoodCorrelationStencil(
  const T *inPtr1, const T *inPtr2,
  const vtkIdType inInc1[3], const vtkIdType inInc2[3],
  const int extent[6], vtkImageStencilData *stencil,
  int radius, int idY, int idZ, U *workPtr)
{
  int xMin = extent[0];
  int xMax = extent[1];
  int r1 = xMin;
  int r2 = xMax;
  int iter = 0;
  int rval = 1;

  // loop over all stencil extents in the row
  do
    {
    int s1 = ((iter == 0) ? xMin : r2 + 1);
    if (stencil)
      {
      rval = stencil->GetNextExtent(r1, r2, xMin, xMax, idY, idZ, iter);
      }
    int s2 = ((rval == 0) ? xMax : r1 - 1);

    if (s1 != s2 + 1)
      {
      // zero everything outside of the stencil
      int k = s2 - s1 + 1;
      inPtr1 += k*inInc1[0];
      inPtr2 += k*inInc2[0];
      do
        {
        workPtr[0] = 0;
        workPtr[1] = 0;
        workPtr[2] = 0;
        workPtr[3] = 0;
        workPtr[4] = 0;
        workPtr[5] = 0;
        workPtr += 6;
        }
      while (--k);
      }

    if (rval == 0)
      {
      break;
      }

    if (r1 != r2 + 1)
      {
      // apply sliding window filter to stencil extent
      int k = r2 - r1 + 1;
      vtkImageNeighborhoodCorrelationX(
        inPtr1, inPtr2, inInc1[0], inInc2[0], k, radius, workPtr);
      workPtr += k*6;
      inPtr1 += k*inInc1[0];
      inPtr2 += k*inInc2[0];
      }
    }
  while (stencil); // don't loop if no stencil
}

//----------------------------------------------------------------------------
// Apply a 2D filter to image slices (either XY or XZ slices).
// The inPtr parameter must be positioned at the correct slice
template<class T, class U>
void vtkImageNeighborhoodCorrelation2D(
  const T *inPtr1, const T *inPtr2,
  const vtkIdType inInc1[3], const vtkIdType inInc2[3],
  const int extent[6], vtkImageStencilData *stencil,
  int radiusX, int radiusZ, int idY,
  U *workPtr, U *rowPtr)
{
  /*  Sliding window

  Each buffer element is one row of partial sums.

  Legend:
      - uninitialized data
      I input data (headPtr)
      B buffered data from previous iterations
      O output data (workPtr)
      P output data from previous iteration

  Input Output Buffer          Stage       Equation
  Index Index  Contents

  0            O----I------- - Initialize  O = I
  1            O----BI------ - Prime       O = O + I
  2            O----BBI----- - Prime       O = O + I
  3            O----BBBI---- - Prime       O = O + I
  4     0      O----BBBBI--- - Prime       O = O + I
  5     1      PO---BBBBBI-- - Lead-in     O = P + I
  6     2      PPO--BBBBBBI- - Lead-in     O = P + I
  7     3      PPPO-BBBBBBBI - Lead-in     O = P + I
  8     4      PPPPOBBBBBBBB I Lead-in     O = P + I
  9     5      PPPPPOBBBBBBB I Slide       O = P - O + I
               PPPPPPOBBBBBB I Slide       O = P - O + I
               PPPPPPPOBBBBB I Slide       O = P - O + I
  n-1   n-5    PPPPPPPPOBBBB I Slide       O = P - O + I
        n-4    PPPPPPPPPOBBB I Finish      O = P - O
        n-3    PPPPPPPPPPOBB I Finish      O = P - O
        n-2    PPPPPPPPPPPOB I Finish      O = P - O
        n-1    PPPPPPPPPPPPO I Finish      O = P - O
  */

  int idZMin = extent[4];
  int idZMax = extent[5];
  int rowSize = extent[1] - extent[0] + 1;
  vtkIdType elementSize = 6; // x,y,xx,yy,xy,n

  U *lastWorkPtr = workPtr;
  U *headPtr = workPtr + elementSize*rowSize*(radiusZ + 1);
  U *endPtr = workPtr + elementSize*rowSize*(idZMax - idZMin + 1);
  U *checkPtr = rowPtr + elementSize*rowSize;

  if (radiusZ + 1 > idZMax - idZMin)
    {
    // use the row buffer
    headPtr = rowPtr;
    }

  // filter in both X and Z directions
  for (int idZ = idZMin; idZ <= idZMax; idZ++)
    {
    if (idZMin == idZMax)
      {
      // if only one row, set headPtr to that row
      headPtr = workPtr;
      }
    else if (headPtr == endPtr || headPtr == checkPtr)
      {
      // use the row buffer when beyond the end of the main buffer
      headPtr = rowPtr;
      }

    // apply the filter in the X direction
    vtkImageNeighborhoodCorrelationStencil(
      inPtr1, inPtr2, inInc1, inInc2, extent, stencil,
      radiusX, idY, idZ, headPtr);

    if (idZMin == idZMax)
      {
      // only one row needed
      return;
      }

    inPtr1 += inInc1[2];
    inPtr2 += inInc2[2];

    if (idZ == idZMin)
      {
      // initialize first row
      U *tmpPtr = workPtr;
      int k = rowSize;
      do
        {
        workPtr[0] = headPtr[0];
        workPtr[1] = headPtr[1];
        workPtr[2] = headPtr[2];
        workPtr[3] = headPtr[3];
        workPtr[4] = headPtr[4];
        workPtr[5] = headPtr[5];
        workPtr += 6;
        headPtr += 6;
        }
      while (--k);
      workPtr = tmpPtr;
      }
    else if (idZ - idZMin < radiusZ + 1)
      {
      // prime the first row
      U *tmpPtr = workPtr;
      int k = rowSize;
      do
        {
        workPtr[0] += headPtr[0];
        workPtr[1] += headPtr[1];
        workPtr[2] += headPtr[2];
        workPtr[3] += headPtr[3];
        workPtr[4] += headPtr[4];
        workPtr[5] += headPtr[5];
        workPtr += 6;
        headPtr += 6;
        }
      while (--k);
      if (idZ - idZMin < radiusZ)
        {
        workPtr = tmpPtr;
        }
      }
    else if (idZ - idZMin < 2*radiusZ + 1)
      {
      // perform the lead-in
      int k = rowSize;
      do
        {
        workPtr[0] = lastWorkPtr[0] + headPtr[0];
        workPtr[1] = lastWorkPtr[1] + headPtr[1];
        workPtr[2] = lastWorkPtr[2] + headPtr[2];
        workPtr[3] = lastWorkPtr[3] + headPtr[3];
        workPtr[4] = lastWorkPtr[4] + headPtr[4];
        workPtr[5] = lastWorkPtr[5] + headPtr[5];
        workPtr += 6;
        lastWorkPtr += 6;
        headPtr += 6;
        }
      while (--k);
      }
    else
      {
      // apply the sliding window
      int k = rowSize;
      do
        {
        workPtr[0] = lastWorkPtr[0] + headPtr[0] - workPtr[0];
        workPtr[1] = lastWorkPtr[1] + headPtr[1] - workPtr[1];
        workPtr[2] = lastWorkPtr[2] + headPtr[2] - workPtr[2];
        workPtr[3] = lastWorkPtr[3] + headPtr[3] - workPtr[3];
        workPtr[4] = lastWorkPtr[4] + headPtr[4] - workPtr[4];
        workPtr[5] = lastWorkPtr[5] + headPtr[5] - workPtr[5];
        workPtr += 6;
        lastWorkPtr += 6;
        headPtr += 6;
        }
      while (--k);
      }
    }

  // if workPtr was never incremented, we're done
  if (workPtr == lastWorkPtr)
    {
    return;
    }

  // finish up the final bit
  for (int i = 0; i < radiusZ; i++)
    {
    if (idZMax - idZMin + i < 2*radiusZ)
      {
      // finishing the "lead-in"
      int k = rowSize;
      do
        {
        workPtr[0] = lastWorkPtr[0];
        workPtr[1] = lastWorkPtr[1];
        workPtr[2] = lastWorkPtr[2];
        workPtr[3] = lastWorkPtr[3];
        workPtr[4] = lastWorkPtr[4];
        workPtr[5] = lastWorkPtr[5];
        workPtr += 6;
        lastWorkPtr += 6;
        }
      while (--k);
      }
    else
      {
      // finishing the "slide"
      int k = rowSize;
      do
        {
        workPtr[0] = lastWorkPtr[0] - workPtr[0];
        workPtr[1] = lastWorkPtr[1] - workPtr[1];
        workPtr[2] = lastWorkPtr[2] - workPtr[2];
        workPtr[3] = lastWorkPtr[3] - workPtr[3];
        workPtr[4] = lastWorkPtr[4] - workPtr[4];
        workPtr[5] = lastWorkPtr[5] - workPtr[5];
        workPtr += 6;
        lastWorkPtr += 6;
        }
      while (--k);
      }
    }
}

//----------------------------------------------------------------------------
// Apply a sliding window filter in all three directions and incrementally
// compute the normalized cross correlation.
template<class T, class U>
void vtkImageNeighborhoodCorrelation3D(
  const T *inPtr1, const T *inPtr2,
  const vtkIdType inInc1[3], const vtkIdType inInc2[3],
  const int extent[6], const int pieceExtent[6],
  vtkImageStencilData *stencil,
  const int radius[3], U *workPtr, vtkAlgorithm *progress,
  vtkImageNeighborhoodCorrelationThreadData *threadLocal)
{
  // apply filter in all three directions: first X, then Z, then Y
  // (doing Z second is most efficient, memory-wise, because it is
  // the dimension broken up between threads)

  /*  Sliding window with rotating buffer:

  Every iteration, the buffer rotates left.
  Each buffer element is a XZ slice of partial sums.

  Legend:
      - uninitialized data
      I input data
      B buffered data from previous iterations
      O output data
      P output data from previous iteration

  Input Output Buffer      Stage       Equation
  Index Index  Contents

  0            -----O----I Initialize  O = I
  1            ----O----BI Prime       O = O + I
  2            ---O----BBI Prime       O = O + I
  3            --O----BBBI Prime       O = O + I
  4     0      -O----BBBBI Prime       O = O + I
  5     1      PO---BBBBBI Lead-in     O = P + I
  6     2      PO--BBBBBBI Lead-in     O = P + I
  7     3      PO-BBBBBBBI Lead-in     O = P + I
  8     4      POBBBBBBBBI Lead-in     O = P + I
  9     5      POBBBBBBBBI Slide       O = P - O + I
               POBBBBBBBBI Slide       O = P - O + I
               POBBBBBBBBI Slide       O = P - O + I
  n-1   n-5    POBBBBBBBBI Slide       O = P - O + I
        n-4    POBBBBBBBB- Finish      O = P - O
        n-3    POBBBBBBB-- Finish      O = P - O
        n-2    POBBBBBB--- Finish      O = P - O
        n-1    POBBBBB---- Finish      O = P - O
  */

  double result = 0.0;

  int radiusX = radius[0];
  int radiusY = radius[1];
  int radiusZ = radius[2];
  int bufferSize = 2*radiusY + 3;

  int idYMin = extent[2];
  int idYMax = extent[3];

  // compute temporary workspace requirements
  vtkIdType elementSize = 6; // x,y,xx,yy,xy,n
  vtkIdType rowSize = extent[1] - extent[0] + 1;
  vtkIdType sliceSize = rowSize*(extent[5] - extent[4] + 1);
  vtkIdType workSize = sliceSize*bufferSize;
  if (extent[5] > extent[4])
    {
    workSize += rowSize;
    }

  // temporary workspace for slices of sums of x,y,xx,yy,xy,n
  U *workPtr2 = new U[workSize*elementSize];
  U **bufferPtr = new U *[bufferSize];
  for (int jj = 0; jj < bufferSize; jj++)
    {
    bufferPtr[jj] = workPtr2 + jj*sliceSize*elementSize;
    }
  // temporary workspace for a single row of x,y,xx,yy,xy,n
  U *rowPtr = workPtr2;
  if (extent[5] > extent[4])
    {
    rowPtr += bufferSize*sliceSize*elementSize;
    }

  // progress reporting variables
  int progressGoal = idYMax - idYMin + 1;
  int progressStep = (progressGoal + 49)/50;
  int progressCount = 0;

  // loop through the XZ slices
  for (int idY = extent[2]; idY <= extent[3]; idY++)
    {
    if (progress != NULL && (progressCount % progressStep) == 0)
      {
      progress->UpdateProgress(progressCount*1.0/progressGoal);
      }
    progressCount++;

    // always put the new data into the first buffer, which
    // will become the last buffer after the buffers are rotated
    U *headPtr = bufferPtr[0];

    // compute the next slice
    vtkImageNeighborhoodCorrelation2D(
      inPtr1, inPtr2, inInc1, inInc2, extent, stencil,
      radiusX, radiusZ, idY, headPtr, rowPtr);

    inPtr1 += inInc1[1];
    inPtr2 += inInc2[1];

    // rotate the buffers
    for (int j = 1; j < bufferSize; j++)
      {
      bufferPtr[j-1] = bufferPtr[j];
      }
    bufferPtr[bufferSize-1] = headPtr;

    if (idY == idYMin)
      {
      // initialize first slice
      workPtr = bufferPtr[radiusY + 1];
      vtkIdType k = sliceSize;
      do
        {
        workPtr[0] = headPtr[0];
        workPtr[1] = headPtr[1];
        workPtr[2] = headPtr[2];
        workPtr[3] = headPtr[3];
        workPtr[4] = headPtr[4];
        workPtr[5] = headPtr[5];
        workPtr += 6;
        headPtr += 6;
        }
      while (--k);
      }
    else if (idY - idYMin < radiusY + 1)
      {
      // prime the first slice
      workPtr = bufferPtr[radiusY + 1 - (idY - idYMin)];
      vtkIdType k = sliceSize;
      do
        {
        workPtr[0] += headPtr[0];
        workPtr[1] += headPtr[1];
        workPtr[2] += headPtr[2];
        workPtr[3] += headPtr[3];
        workPtr[4] += headPtr[4];
        workPtr[5] += headPtr[5];
        workPtr += 6;
        headPtr += 6;
        }
      while (--k);
      }
    else if (idY - idYMin < 2*radiusY + 1)
      {
      // perform the lead-in
      U *lastWorkPtr = bufferPtr[0];
      workPtr = bufferPtr[1];
      vtkIdType k = sliceSize;
      do
        {
        workPtr[0] = lastWorkPtr[0] + headPtr[0];
        workPtr[1] = lastWorkPtr[1] + headPtr[1];
        workPtr[2] = lastWorkPtr[2] + headPtr[2];
        workPtr[3] = lastWorkPtr[3] + headPtr[3];
        workPtr[4] = lastWorkPtr[4] + headPtr[4];
        workPtr[5] = lastWorkPtr[5] + headPtr[5];
        workPtr += 6;
        lastWorkPtr += 6;
        headPtr += 6;
        }
      while (--k);
      }
    else
      {
      // apply the sliding window
      U *lastWorkPtr = bufferPtr[0];
      workPtr = bufferPtr[1];
      vtkIdType k = sliceSize;
      do
        {
        workPtr[0] = lastWorkPtr[0] + headPtr[0] - workPtr[0];
        workPtr[1] = lastWorkPtr[1] + headPtr[1] - workPtr[1];
        workPtr[2] = lastWorkPtr[2] + headPtr[2] - workPtr[2];
        workPtr[3] = lastWorkPtr[3] + headPtr[3] - workPtr[3];
        workPtr[4] = lastWorkPtr[4] + headPtr[4] - workPtr[4];
        workPtr[5] = lastWorkPtr[5] + headPtr[5] - workPtr[5];
        workPtr += 6;
        lastWorkPtr += 6;
        headPtr += 6;
        }
      while (--k);
      }

    // read as "if (idY - idYMin >= radiusY) { for (;;) { ... } }",
    // and see break statement below.
    int i = 0;
    while (idY - idYMin >= radiusY)
      {
      // the sums over the neighborhoods have been computed for all the
      // voxels in a slice, so compute the normalized cross-correlation
      // (only compute the metric over the pieceExtent)
      int outIdY = idY - radiusY + i;
      if (outIdY >= pieceExtent[2] && outIdY <= pieceExtent[3])
        {
        double total = 0;
        workPtr = bufferPtr[1];
        workPtr += elementSize*rowSize*(pieceExtent[4] - extent[4]);
        for (int idZ = pieceExtent[4]; idZ <= pieceExtent[5]; idZ++)
          {
          workPtr += elementSize*(pieceExtent[0] - extent[0]);

          // only compute the metric within the stencil
          int iter = 0;
          int rval = 1;
          int r1 = pieceExtent[0];
          int r2 = pieceExtent[1];

          // loop over stencil extents (break at end if no stencil)
          do
            {
            int s1 = ((iter == 0) ? pieceExtent[0] : r2 + 1);
            if (stencil)
              {
              rval = stencil->GetNextExtent(
                r1, r2, pieceExtent[0], pieceExtent[1], outIdY, idZ, iter);
              }
            int s2 = ((rval == 0) ? pieceExtent[1] : r1 - 1);
            workPtr += elementSize*(s2 - s1 + 1);

            if (rval == 0)
              {
              break;
              }

            if (r1 != r2 + 1)
              {
              int kk = r2 - r1 + 1;
              do
                {
                U xSum = workPtr[0];
                U ySum = workPtr[1];
                U xxSum = workPtr[2];
                U yySum = workPtr[3];
                U xySum = workPtr[4];
                U count = workPtr[5];
                workPtr += 6;
                double denom = static_cast<double>(xxSum*count - xSum*xSum)*
                  static_cast<double>(yySum*count - ySum*ySum);
                double numer = static_cast<double>(xySum*count - xSum*ySum);
                numer *= numer;
                double nccSquared = 1.0;
                if (denom > 0)
                  {
                  nccSquared = numer/denom;
                  }
                total += nccSquared;
                }
              while (--kk);
              }
            }
          while (stencil);

          workPtr += elementSize*(extent[1] - pieceExtent[1]);
          }

        workPtr += elementSize*rowSize*(extent[5] - pieceExtent[5]);

        result += total;
        }

      if (idY < idYMax || outIdY >= pieceExtent[3])
        {
        break;
        }

      // rotate the buffers
      headPtr = bufferPtr[0];
      for (int j = 1; j < bufferSize; j++)
        {
        bufferPtr[j-1] = bufferPtr[j];
        }
      bufferPtr[bufferSize-1] = headPtr;

      // tail off if there is any tail left
      U *lastWorkPtr = bufferPtr[0];
      workPtr = bufferPtr[1];
      if (idYMax - idYMin + i < 2*radiusY)
        {
        // finish lead-in
        vtkIdType k = sliceSize;
        do
          {
          workPtr[0] = lastWorkPtr[0];
          workPtr[1] = lastWorkPtr[1];
          workPtr[2] = lastWorkPtr[2];
          workPtr[3] = lastWorkPtr[3];
          workPtr[4] = lastWorkPtr[4];
          workPtr[5] = lastWorkPtr[5];
          workPtr += 6;
          lastWorkPtr += 6;
          }
        while (--k);
        }
      else
        {
        // finish slide
        vtkIdType k = sliceSize;
        do
          {
          workPtr[0] = lastWorkPtr[0] - workPtr[0];
          workPtr[1] = lastWorkPtr[1] - workPtr[1];
          workPtr[2] = lastWorkPtr[2] - workPtr[2];
          workPtr[3] = lastWorkPtr[3] - workPtr[3];
          workPtr[4] = lastWorkPtr[4] - workPtr[4];
          workPtr[5] = lastWorkPtr[5] - workPtr[5];
          workPtr += 6;
          lastWorkPtr += 6;
          }
        while (--k);
        }

      // counter for tail loop
      i++;
      }
    }

  delete [] workPtr2;
  delete [] bufferPtr;

  threadLocal->Result += result;
}

} // end anonymous namespace

//----------------------------------------------------------------------------
int vtkImageNeighborhoodCorrelation::RequestData(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // create the thread-local object
  vtkImageNeighborhoodCorrelationTLS tlocal;
  tlocal.Initialize(this);
  this->ThreadData = &tlocal;

  this->Superclass::RequestData(request, inputVector, outputVector);

  this->ThreadData = 0;

  return 1;
}
//----------------------------------------------------------------------------
// This method is passed a input and output region, and executes the filter
// algorithm to fill the output from the input.
// It just executes a switch statement to call the correct function for
// the regions data types.
void vtkImageNeighborhoodCorrelation::PieceRequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector),
  const int pieceExtent[6], vtkIdType pieceId)
{
  vtkInformation *inInfo0 = inputVector[0]->GetInformationObject(0);
  vtkInformation *inInfo1 = inputVector[1]->GetInformationObject(0);

  vtkImageData *inData0 = vtkImageData::SafeDownCast(
    inInfo0->Get(vtkDataObject::DATA_OBJECT()));
  vtkImageData *inData1 = vtkImageData::SafeDownCast(
    inInfo1->Get(vtkDataObject::DATA_OBJECT()));

  if (inData0->GetScalarType() != inData1->GetScalarType())
    {
    if (pieceId == 0)
      {
      vtkErrorMacro("input image types must be the same.");
      }
    return;
    }

  // copy this, in case the user changes it during execution
  int neighborhoodRadius[3];
  neighborhoodRadius[0] = this->NeighborhoodRadius[0];
  neighborhoodRadius[1] = this->NeighborhoodRadius[1];
  neighborhoodRadius[2] = this->NeighborhoodRadius[2];

  // make sure execute extent is not beyond the extent of any input
  int inExt0[6], inExt1[6];
  inData0->GetExtent(inExt0);
  inData1->GetExtent(inExt1);

  int extent[6], pExtent[6];
  for (int i = 0; i < 6; i += 2)
    {
    // partial sums need to be computed over pieceExtent + radius, to
    // ensure there is no "gap" between slabs assigned to different threads
    int j = i + 1;
    int r = neighborhoodRadius[i/2];
    extent[i] = pieceExtent[i] - r;
    extent[j] = pieceExtent[j] + r;
    extent[i] = ((extent[i] > inExt0[i]) ? extent[i] : inExt0[i]);
    extent[i] = ((extent[i] > inExt1[i]) ? extent[i] : inExt1[i]);
    extent[j] = ((extent[j] < inExt0[j]) ? extent[j] : inExt0[j]);
    extent[j] = ((extent[j] < inExt1[j]) ? extent[j] : inExt1[j]);
    pExtent[i] =
      ((pieceExtent[i] > extent[i]) ? pieceExtent[i] : extent[i]);
    pExtent[j] =
      ((pieceExtent[j] < extent[j]) ? pieceExtent[j] : extent[j]);
    if (pExtent[i] > pExtent[j])
      {
      return;
      }
    }

  void *inPtr0 = inData0->GetScalarPointerForExtent(extent);
  void *inPtr1 = inData1->GetScalarPointerForExtent(extent);

  vtkIdType inInc1[3], inInc2[3];
  inData0->GetIncrements(inInc1);
  inData1->GetIncrements(inInc2);

  vtkImageStencilData *stencil = this->GetStencil();

  // only used for tracking progress
  vtkAlgorithm *progress = (pieceId == 0 ? this : 0);

  int scalarType = inData0->GetScalarType();

  if (scalarType == VTK_FLOAT || scalarType == VTK_DOUBLE)
    {
    // use a floating-point type for computing sums
    double workVal = 0;

    if (scalarType == VTK_FLOAT)
      {
      vtkImageNeighborhoodCorrelation3D(
        static_cast<float *>(inPtr0), static_cast<float *>(inPtr1),
        inInc1, inInc2, extent, pExtent, stencil, neighborhoodRadius,
        &workVal, progress, &this->ThreadData->Local(pieceId));
      }
    else
      {
      vtkImageNeighborhoodCorrelation3D(
        static_cast<double *>(inPtr0), static_cast<double *>(inPtr1),
        inInc1, inInc2, extent, pExtent, stencil, neighborhoodRadius,
        &workVal, progress, &this->ThreadData->Local(pieceId));
      }
    }
  else
    {
    // use an integer type for computing sums
    vtkTypeInt64 workVal = 0;

    // turn off floats in the vtkTemplateAliasMacro
#undef VTK_USE_FLOAT64
#define VTK_USE_FLOAT64 0
#undef VTK_USE_FLOAT32
#define VTK_USE_FLOAT32 0

    switch (scalarType)
      {
      vtkTemplateAliasMacro(
        vtkImageNeighborhoodCorrelation3D(
          static_cast<VTK_TT *>(inPtr0), static_cast<VTK_TT *>(inPtr1),
          inInc1, inInc2, extent, pExtent, stencil, neighborhoodRadius,
          &workVal, progress, &this->ThreadData->Local(pieceId)));
      default:
        vtkErrorMacro(<< "Execute: Unknown ScalarType");
      }
    }
}

//----------------------------------------------------------------------------
void vtkImageNeighborhoodCorrelation::ReduceRequestData(
  vtkInformation *, vtkInformationVector **, vtkInformationVector *)
{
  double result = 0.0;

  for (vtkImageNeighborhoodCorrelationTLS::iterator
       iter = this->ThreadData->begin();
       iter != this->ThreadData->end(); ++iter)
    {
    result += iter->Result;
    }

  this->SetMinimizable(-result);
}
