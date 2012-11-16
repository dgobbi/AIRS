/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageIslandRemoval.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkImageIslandRemoval.h"

#include "vtkMath.h"
#include "vtkImageData.h"
#include "vtkImageStencilData.h"
#include "vtkImageRegionIterator.h"
#include "vtkImageIterator.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkImageStencilData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkTemplateAliasMacro.h"

#include <stack>
#include <algorithm>

vtkStandardNewMacro(vtkImageIslandRemoval);

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageIslandRemoval::vtkImageIslandRemoval()
{
  this->UpperThreshold = VTK_LARGE_FLOAT;
  this->LowerThreshold = -VTK_LARGE_FLOAT;
  this->ReplaceIn = 0;
  this->InValue = 0.0;
  this->ReplaceOut = 0;
  this->OutValue = 0.0;

  this->NeighborhoodRadius[0] = 0.0;
  this->NeighborhoodRadius[1] = 0.0;
  this->NeighborhoodRadius[2] = 0.0;
  this->NeighborhoodFraction = 0.5;

  this->SliceRangeX[0] = -VTK_INT_MAX;
  this->SliceRangeX[1] = VTK_INT_MAX;
  this->SliceRangeY[0] = -VTK_INT_MAX;
  this->SliceRangeY[1] = VTK_INT_MAX;
  this->SliceRangeZ[0] = -VTK_INT_MAX;
  this->SliceRangeZ[1] = VTK_INT_MAX;

  this->ActiveComponent = -1;

  this->ImageMask = vtkImageData::New();

  this->NumberOfInVoxels = 0;

  this->SetNumberOfInputPorts(2);
}

//----------------------------------------------------------------------------
vtkImageIslandRemoval::~vtkImageIslandRemoval()
{
  this->ImageMask->Delete();
}

//----------------------------------------------------------------------------
void vtkImageIslandRemoval::SetInValue(double val)
{
  if (val != this->InValue || this->ReplaceIn != 1)
    {
    this->InValue = val;
    this->ReplaceIn = 1;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkImageIslandRemoval::SetOutValue(double val)
{
  if (val != this->OutValue || this->ReplaceOut != 1)
    {
    this->OutValue = val;
    this->ReplaceOut = 1;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
// The values greater than or equal to the value match.
void vtkImageIslandRemoval::ThresholdByUpper(double thresh)
{
  if (this->LowerThreshold != thresh || this->UpperThreshold < VTK_LARGE_FLOAT)
    {
    this->LowerThreshold = thresh;
    this->UpperThreshold = VTK_LARGE_FLOAT;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
// The values less than or equal to the value match.
void vtkImageIslandRemoval::ThresholdByLower(
  double thresh)
{
  if (this->UpperThreshold != thresh ||
      this->LowerThreshold > -VTK_LARGE_FLOAT)
    {
    this->UpperThreshold = thresh;
    this->LowerThreshold = -VTK_LARGE_FLOAT;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
// The values in a range (inclusive) match
void vtkImageIslandRemoval::ThresholdBetween(
  double lower, double upper)
{
  if (this->LowerThreshold != lower ||
      this->UpperThreshold != upper)
    {
    this->LowerThreshold = lower;
    this->UpperThreshold = upper;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
int vtkImageIslandRemoval::FillInputPortInformation(
  int port, vtkInformation* info)
{
  if (port == 1)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageStencilData");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }
  else
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    }
  return 1;
}

//----------------------------------------------------------------------------
void vtkImageIslandRemoval::SetStencil(vtkImageStencilData *stencil)
{
  this->SetInput(1, stencil);
}

//----------------------------------------------------------------------------
vtkImageStencilData *vtkImageIslandRemoval::GetStencil()
{
  if (this->GetNumberOfInputConnections(1) < 1)
    {
    return NULL;
    }

  return vtkImageStencilData::SafeDownCast(
    this->GetExecutive()->GetInputData(1, 0));
}

//----------------------------------------------------------------------------
unsigned long vtkImageIslandRemoval::GetMTime()
{
  unsigned long mTime = this->MTime.GetMTime();

  return mTime;
}

//----------------------------------------------------------------------------
// seed struct: just a set of indices
namespace {

class vtkFloodFillSeed
{
public:
  vtkFloodFillSeed() {
    store[0]=0; store[1]=0; store[2]=0; };
  vtkFloodFillSeed(int i, int j, int k) {
    store[0]=i; store[1]=j; store[2]=k; };
  vtkFloodFillSeed(const vtkFloodFillSeed &seed) {
    store[0]=seed.store[0]; store[1]=seed.store[1]; store[2]=seed.store[2]; };
  const int &operator[](int i) const { return store[i]; };
  const vtkFloodFillSeed &operator=(const vtkFloodFillSeed seed) {
    store[0]=seed.store[0]; store[1]=seed.store[1]; store[2]=seed.store[2];
    return *this; };

private:
  int store[3];
};

}

//----------------------------------------------------------------------------
// Make sure the thresholds are valid for the input scalar range
template <class IT>
void vtkImageIslandRemovalThresholds(
  vtkImageIslandRemoval *self,
  vtkImageData *inData, IT &lowerThreshold, IT &upperThreshold)
{
  if (self->GetLowerThreshold() < inData->GetScalarTypeMin())
    {
    lowerThreshold = static_cast<IT>(inData->GetScalarTypeMin());
    }
  else
    {
    if (self->GetLowerThreshold() > inData->GetScalarTypeMax())
      {
      lowerThreshold = static_cast<IT>(inData->GetScalarTypeMax());
      }
    else
      {
      lowerThreshold = static_cast<IT>(self->GetLowerThreshold());
      }
    }
  if (self->GetUpperThreshold() > inData->GetScalarTypeMax())
    {
    upperThreshold = static_cast<IT>(inData->GetScalarTypeMax());
    }
  else
    {
    if (self->GetUpperThreshold() < inData->GetScalarTypeMin())
      {
      upperThreshold = static_cast<IT>(inData->GetScalarTypeMin());
      }
    else
      {
      upperThreshold = static_cast<IT>(self->GetUpperThreshold());
      }
    }
}

//----------------------------------------------------------------------------
// Make sure the replacement values are within the output scalar range
template <class OT>
void vtkImageIslandRemovalValues(
  vtkImageIslandRemoval *self,
  vtkImageData *outData, OT &inValue, OT &outValue)
{
  if (self->GetInValue() < outData->GetScalarTypeMin())
    {
    inValue = static_cast<OT>(outData->GetScalarTypeMin());
    }
  else
    {
    if (self->GetInValue() > outData->GetScalarTypeMax())
      {
      inValue = static_cast<OT>(outData->GetScalarTypeMax());
      }
    else
      {
      inValue = static_cast<OT>(self->GetInValue());
      }
    }
  if (self->GetOutValue() > outData->GetScalarTypeMax())
    {
    outValue = static_cast<OT>(outData->GetScalarTypeMax());
    }
  else
    {
    if (self->GetOutValue() < outData->GetScalarTypeMin())
      {
      outValue = static_cast<OT>(outData->GetScalarTypeMin());
      }
    else
      {
      outValue = static_cast<OT>(self->GetOutValue());
      }
    }
}

//----------------------------------------------------------------------------
void vtkImageIslandRemovalApplyStencil(
  vtkImageData *maskData, vtkImageStencilData *stencil, int extent[6])
{
  vtkImageRegionIterator<unsigned char> iter(maskData, stencil, extent);
  while (!iter.IsAtEnd())
    {
    unsigned char *beginptr = iter.BeginSpan();
    unsigned char *endptr = iter.EndSpan();
    unsigned char val = (iter.IsInStencil() ? 0 : 1);

    for (unsigned char *ptr = beginptr; ptr < endptr; ptr++)
      {
      *ptr = val;
      }

    iter.NextSpan();
    }
}

//----------------------------------------------------------------------------
struct vtkImageIslandRemovalCompare
{
  vtkIdType *Array;

  vtkImageIslandRemovalCompare(vtkIdType *a) : Array(a) {};

  bool operator()(unsigned char x, unsigned char y)
  {
    return (Array[x] > Array[y]);
  }
};

//----------------------------------------------------------------------------
void vtkImageIslandRemovalSort(
  vtkIdType islandSize[256], unsigned char indexed[256])
{
  unsigned char sorted[256];
  vtkIdType oldSize[256];
  for (int i = 0; i < 256; i++)
    {
    sorted[i] = i;
    oldSize[i] = islandSize[i];
    }

  vtkImageIslandRemovalCompare cmpfunc(oldSize);
  std::sort(&sorted[1], &sorted[256], cmpfunc);

  for (int i = 0; i < 256; i++)
    {
    indexed[sorted[i]] = i;
    islandSize[i] = oldSize[sorted[i]];
    }
}

//----------------------------------------------------------------------------
// prune all but the largest "n" islands
void vtkImageIslandRemovalPrune(
  vtkImageData *maskData, vtkImageStencilData *stencil, int extent[6],
  vtkIdType islandSize[256], int n)
{
  unsigned char u = static_cast<unsigned char>(n + 1);

  unsigned char indexed[256];
  vtkImageIslandRemovalSort(islandSize, indexed);

  vtkImageRegionIterator<unsigned char> iter(maskData, stencil, extent);

  while (!iter.IsAtEnd())
    {
    unsigned char *beginptr = iter.BeginSpan();
    unsigned char *endptr = iter.EndSpan();

    if (iter.IsInStencil())
      {
      for (unsigned char *ptr = beginptr; ptr != endptr; ptr++)
        {
        unsigned char v = indexed[*ptr];
        *ptr = (v < u ? v : u);
        }
      }

    iter.NextSpan();
    }
}

//----------------------------------------------------------------------------
template<class IT, class OT>
void vtkImageIslandRemovalFinish(vtkImageIslandRemoval *self,
  vtkImageData *maskData, vtkImageStencilData *stencil, int extent[6],
  vtkImageData *inData, IT *inPtr, vtkImageData *outData, OT *outPtr,
  vtkIdType islandSize[256], int n)
{
  unsigned char u = static_cast<unsigned char>(n + 1);

  unsigned char indexed[256];
  vtkImageIslandRemovalSort(islandSize, indexed);

  // Get replace values as output data type
  bool replaceIn = (self->GetReplaceIn() != 0);
  bool replaceOut = (self->GetReplaceOut() != 0);
  OT inValue, outValue;
  vtkImageIslandRemovalValues(self, outData, inValue, outValue);

  vtkImageIterator<IT> inIt(inData, extent);
  vtkImageIterator<OT> outIt(outData, extent);
  vtkImageRegionIterator<unsigned char> iter(maskData, stencil, extent);

  inPtr = inIt.BeginSpan();
  outPtr = outIt.BeginSpan();
  while (!iter.IsAtEnd())
    {
    unsigned char *ptr = iter.BeginSpan();
    size_t count = iter.EndSpan() - ptr;

    if (iter.IsInStencil())
      {
      if (replaceIn && replaceOut)
        {
        while (count)
          {
          unsigned char v = indexed[*ptr++];
          v = (v < u ? v : 0);
          *outPtr++ = (v != 0 ? inValue : outValue);
          --count;
          }
        }
      else if (!replaceIn && replaceOut)
        {
        while (count)
          {
          unsigned char v = indexed[*ptr++];
          v = (v < u ? v : 0);
          *outPtr++ = (v != 0 ? *inPtr : outValue);
          inPtr++;
          --count;
          }
        }
      else if (replaceIn && !replaceOut)
        {
        while (count)
          {
          unsigned char v = indexed[*ptr++];
          v = (v < u ? v : 0);
          *outPtr++ = (v != 0 ? inValue : *inPtr);
          inPtr++;
          --count;
          }
        }
      else
        {
        while (count)
          {
          *outPtr++ = *inPtr++;
          --count;
          }
        }
      }
    else
      {
      if (replaceOut)
        {
        while (count)
          {
          *outPtr++ = outValue;
          --count;
          }
        }
      else
        {
        while (count)
          {
          *outPtr++ = *inPtr++;
          --count;
          }
        }
      }

    iter.NextSpan();
    if (outPtr == outIt.EndSpan())
      {
      outIt.NextSpan();
      outPtr = outIt.BeginSpan();
      inIt.NextSpan();
      inPtr = inIt.BeginSpan();
      }
    }
}

//----------------------------------------------------------------------------
// This templated function executes the filter for any type of data.
template <class IT, class OT>
void vtkImageIslandRemovalExecute(
  vtkImageIslandRemoval *self,
  vtkImageData *inData, vtkImageData *outData, vtkImageStencilData *stencil,
  vtkImageData *maskData, int outExt[6], int id, IT *inPtr, OT *outPtr,
  int &voxelCount)
{
  // Get active component (only one component is thresholded)
  int nComponents = outData->GetNumberOfScalarComponents();
  int activeComponent = self->GetActiveComponent();
  if (activeComponent < 0) { activeComponent = 0; }
  activeComponent = activeComponent % nComponents;

  // Get thresholds as input data type
  IT lowerThreshold, upperThreshold;
  vtkImageIslandRemovalThresholds(
    self, inData, lowerThreshold, upperThreshold);

  // Get the extent for the flood fill, and clip with the input extent
  int extent[6];
  int inExt[6];
  self->GetSliceRangeX(&extent[0]);
  self->GetSliceRangeY(&extent[2]);
  self->GetSliceRangeZ(&extent[4]);
  inData->GetExtent(inExt);
  int outCheck = 0;
  for (int ii = 0; ii < 3; ii++)
    {
    if (extent[2*ii] > inExt[2*ii+1] || extent[2*ii+1] < inExt[2*ii])
      { // extents don't intersect, we're done
      return;
      }
    if (extent[2*ii] < inExt[2*ii])
      {
      extent[2*ii] = inExt[2*ii];
      }
    if (extent[2*ii+1] > inExt[2*ii+1])
      {
      extent[2*ii+1] = inExt[2*ii+1];
      }
    // check against output extent
    if (extent[2*ii] < outExt[2*ii] || extent[2*ii+1] > outExt[2*ii+1])
      {
      outCheck = 1;
      }
    }

  // Indexing goes from 0 to maxIdX
  int maxIdX = extent[1] - extent[0];
  int maxIdY = extent[3] - extent[2];
  int maxIdZ = extent[5] - extent[4];

  // Total number of voxels
  vtkIdType fullsize = (static_cast<vtkIdType>(extent[1] - extent[0] + 1)*
                        static_cast<vtkIdType>(extent[3] - extent[2] + 1)*
                        static_cast<vtkIdType>(extent[5] - extent[4] + 1));

  // for the progress meter
  double progress = 0.0;
  vtkIdType target = static_cast<vtkIdType>(fullsize/50.0);
  target++;

  // Setup the mask
  maskData->SetOrigin(inData->GetOrigin());
  maskData->SetSpacing(inData->GetSpacing());
  maskData->SetExtent(extent);
  maskData->AllocateScalars(/*VTK_UNSIGNED_CHAR, 1*/);

  unsigned char *maskPtr =
    static_cast<unsigned char *>(maskData->GetScalarPointerForExtent(extent));

  // Get pointers for the new extent (the one used for the seeds)
  inPtr = static_cast<IT *>(inData->GetScalarPointerForExtent(extent));
  outPtr = static_cast<OT *>(outData->GetScalarPointerForExtent(extent));

  // Adjust pointers to active component
  inPtr += activeComponent;
  outPtr += activeComponent;

  if (stencil == 0)
    {
    memset(maskPtr, 0, fullsize);
    }
  else
    {
    // pre-set all mask voxels that are outside the stencil
    vtkImageIslandRemovalApplyStencil(maskData, stencil, extent);
    }

  // Check whether neighborhood will be used
  double f = self->GetNeighborhoodFraction();
  double radius[3];
  self->GetNeighborhoodRadius(radius);
  int xradius = static_cast<int>(radius[0] + 0.5);
  int yradius = static_cast<int>(radius[1] + 0.5);
  int zradius = static_cast<int>(radius[2] + 0.5);
  double fx = 0.0, fy = 0.0, fz = 0.0;
  bool useNeighborhood = ((xradius > 0) & (yradius > 0) & (zradius > 0));
  if (useNeighborhood)
    {
    fx = 1.0/radius[0];
    fy = 1.0/radius[1];
    fz = 1.0/radius[2];
    }

  // Perform the flood fill within the extent
  vtkIdType inInc[3];
  vtkIdType outInc[3];
  vtkIdType maskInc[3];
  inData->GetIncrements(inInc);
  outData->GetIncrements(outInc);
  maskInc[0] = 1;
  maskInc[1] = (extent[1] - extent[0] + 1);
  maskInc[2] = maskInc[1]*(extent[3] - extent[2] + 1);

  double spacing[3];
  double origin[3];
  outData->GetSpacing(spacing);
  outData->GetOrigin(origin);

  // create the seed stack:
  // stack has methods empty(), top(), push(), and pop()
  std::stack<vtkFloodFillSeed> seedStack;

  vtkIdType counter = 0;
  vtkIdType fullcount = 0;

  vtkImageRegionIterator<unsigned char> iter(maskData, stencil, extent);

  unsigned char maskVal = 0;
  vtkIdType islandSize[256];

  while (!iter.IsAtEnd())
    {
    unsigned char *beginptr = iter.BeginSpan();
    unsigned char *endptr = iter.EndSpan();

    if (iter.IsInStencil())
      {
      int xIdx = iter.GetIndexX();
      int yIdx = iter.GetIndexY();
      int zIdx = iter.GetIndexZ();

      IT *inPtr0 = inPtr + (xIdx*inInc[0] + yIdx*inInc[1] + zIdx*inInc[2]);

      for (unsigned char *ptr = beginptr; ptr < endptr; ptr++)
        {
        IT tval = *inPtr0;
        if (*ptr == 0 && lowerThreshold <= tval && tval <= upperThreshold)
          {
          seedStack.push(
            vtkFloodFillSeed(xIdx + extent[0],
                             yIdx + extent[2],
                             zIdx + extent[4]));

          maskVal++;
          if (maskVal == 0)
            {
            // need to prune all but N islands
            vtkImageIslandRemovalPrune(
              maskData, stencil, extent, islandSize, 3);
            maskVal = 5;
            }
          islandSize[maskVal] = 0;

  //====
  while (!seedStack.empty())
    {
    vtkFloodFillSeed seed = seedStack.top();
    seedStack.pop();

    unsigned char *maskPtr1 = maskPtr + (seed[0]*maskInc[0] +
                                         seed[1]*maskInc[1] +
                                         seed[2]*maskInc[2]);

    if (*maskPtr1)
      {
      continue;
      }
    *maskPtr1 = maskVal;

    fullcount++;
    if (id == 0 && (fullcount % target) == 0)
      {
      double v = counter/(10.0*fullcount);
      double p = fullcount/(v*fullsize + (1.0 - v)*fullcount);
      if (p > progress)
        {
        progress = p;
        self->UpdateProgress(progress);
        }
      }

    IT *inPtr1 = inPtr + (seed[0]*inInc[0] +
                          seed[1]*inInc[1] +
                          seed[2]*inInc[2]);
    IT temp = *inPtr1;

    bool inside = ((lowerThreshold <= temp) & (temp <= upperThreshold));

    // use a spherical neighborhood
    if (useNeighborhood)
      {
      int xmin = seed[0] - xradius;
      xmin = (xmin >= 0 ? xmin : 0);
      int xmax = seed[0] + xradius;
      xmax = (xmax <= maxIdX ? xmax : maxIdX);

      int ymin = seed[1] - yradius;
      ymin = (ymin >= 0 ? ymin : 0);
      int ymax = seed[1] + yradius;
      ymax = (ymax <= maxIdY ? ymax : maxIdY);

      int zmin = seed[2] - zradius;
      zmin = (zmin >= 0 ? zmin : 0);
      int zmax = seed[2] + zradius;
      zmax = (zmax <= maxIdZ ? zmax : maxIdZ);

      inPtr1 = inPtr + (xmin*inInc[0] +
                        ymin*inInc[1] +
                        zmin*inInc[2]);

      int totalcount = 0;
      int threshcount = 0;
      int iz = zmin;
      do
        {
        IT *inPtr2 = inPtr1;
        double rz = (iz - seed[2])*fz;
        rz *= rz;
        int iy = ymin;
        do
          {
          IT *inPtr3 = inPtr2;
          double ry = (iy - seed[1])*fy;
          ry *= ry;
          double rzy = rz + ry;
          int ix = xmin;
          do
            {
            double rx = (ix - seed[0])*fx;
            rx *= rx;
            double rzyx = rzy + rx;
            // include a tolerance in radius check
            bool isgood = (rzyx < (1.0 + 7.62939453125e-06));
            totalcount += isgood;
            isgood &= ((lowerThreshold <= *inPtr3) &
                       (*inPtr3 <= upperThreshold));
            threshcount += isgood;
            inPtr3 += inInc[0];
            }
          while (++ix <= xmax);
          inPtr2 += inInc[1];
          }
        while (++iy <= ymax);
        inPtr1 += inInc[2];
        }
      while (++iz <= zmax);

      // what fraction of the sphere is within threshold?
      inside &= !(static_cast<double>(threshcount) < totalcount*f);
      }

    if (inside)
      {
      // count the seed
      counter += 1;
      islandSize[maskVal] += 1;

      // push the new seeds
      if (seed[2] > 0 && *(maskPtr1 - maskInc[2]) == 0)
        {
        seedStack.push(vtkFloodFillSeed(seed[0], seed[1], seed[2] - 1));
        }
      if (seed[2] < maxIdZ && *(maskPtr1 + maskInc[2]) == 0)
        {
        seedStack.push(vtkFloodFillSeed(seed[0], seed[1], seed[2] + 1));
        }
      if (seed[1] > 0 && *(maskPtr1 - maskInc[1]) == 0)
        {
        seedStack.push(vtkFloodFillSeed(seed[0], seed[1] - 1, seed[2]));
        }
      if (seed[1] < maxIdY && *(maskPtr1 + maskInc[1]) == 0)
        {
        seedStack.push(vtkFloodFillSeed(seed[0], seed[1] + 1, seed[2]));
        }
      if (seed[0] > 0 && *(maskPtr1 - maskInc[0]) == 0)
        {
        seedStack.push(vtkFloodFillSeed(seed[0] - 1, seed[1], seed[2]));
        }
      if (seed[0] < maxIdX && *(maskPtr1 + maskInc[0]) == 0)
        {
        seedStack.push(vtkFloodFillSeed(seed[0] + 1, seed[1], seed[2]));
        }
      }
    }
  //====

          }

        inPtr0 += inInc[0];
        xIdx++;
        }
      }
    iter.NextSpan();
    }

  vtkImageIslandRemovalFinish(
    self, maskData, stencil, extent,
    inData, inPtr, outData, outPtr, islandSize, 1);

  if (id == 0)
    {
    self->UpdateProgress(1.0);
    }

  voxelCount = counter;
}

//----------------------------------------------------------------------------
int vtkImageIslandRemoval::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector))
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *stencilInfo = inputVector[1]->GetInformationObject(0);

  int inExt[6], extent[6];
  inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), inExt);
  extent[0] = this->SliceRangeX[0];
  extent[1] = this->SliceRangeX[1];
  extent[2] = this->SliceRangeY[0];
  extent[3] = this->SliceRangeY[1];
  extent[4] = this->SliceRangeZ[0];
  extent[5] = this->SliceRangeZ[1];

  // Clip the extent to the inExt
  for (int i = 0; i < 3; i++)
    {
    if (extent[2*i] < inExt[2*i])
      {
      extent[2*i] = inExt[2*i];
      }
    if (extent[2*i+1] > inExt[2*i+1])
      {
      extent[2*i+1] = inExt[2*i+1];
      }
    }

  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent, 6);
  if (stencilInfo)
    {
    stencilInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
                     extent, 6);
    }

  return 1;
}

//----------------------------------------------------------------------------
int vtkImageIslandRemoval::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *stencilInfo = inputVector[1]->GetInformationObject(0);

  vtkImageData* outData = static_cast<vtkImageData *>(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkImageData* inData = static_cast<vtkImageData *>(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkImageData *maskData = this->ImageMask;

  vtkImageStencilData* stencil = 0;
  if (stencilInfo)
    {
    stencil = static_cast<vtkImageStencilData *>(
      stencilInfo->Get(vtkDataObject::DATA_OBJECT()));
    }

  int outExt[6];
  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), outExt);
  this->AllocateOutputData(outData/*, outInfo*/, outExt);

  // get scalar pointers
  void *inPtr = inData->GetScalarPointerForExtent(outExt);
  void *outPtr = outData->GetScalarPointerForExtent(outExt);

  int id = 0; // not multi-threaded

  if (inData->GetScalarType() != outData->GetScalarType())
    {
    vtkErrorMacro("Execute: Output ScalarType "
                  << outData->GetScalarType()
                  << ", must Input ScalarType "
                  << inData->GetScalarType());
    return 0;
    }

  switch (inData->GetScalarType())
    {
    vtkTemplateAliasMacro(
      vtkImageIslandRemovalExecute(
        this, inData, outData, stencil, maskData, outExt, id,
        static_cast<VTK_TT *>(inPtr), static_cast<VTK_TT *>(outPtr),
        this->NumberOfInVoxels));

    default:
      vtkErrorMacro(<< "Execute: Unknown input ScalarType");
      return 0;
    }

  return 1;
}

//----------------------------------------------------------------------------
void vtkImageIslandRemoval::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "InValue: " << this->InValue << "\n";
  os << indent << "OutValue: " << this->OutValue << "\n";
  os << indent << "LowerThreshold: " << this->LowerThreshold << "\n";
  os << indent << "UpperThreshold: " << this->UpperThreshold << "\n";
  os << indent << "ReplaceIn: " << this->ReplaceIn << "\n";
  os << indent << "ReplaceOut: " << this->ReplaceOut << "\n";
  os << indent << "NeighborhoodRadius: " << this->NeighborhoodRadius[0] << " "
     << this->NeighborhoodRadius[1] << " "
     << this->NeighborhoodRadius[2] << "\n";
  os << indent << "NeighborhoodFraction: "
     << this->NeighborhoodFraction << "\n";
  os << indent << "NumberOfInVoxels: " << this->NumberOfInVoxels << "\n";
  os << indent << "SliceRangeX: "
     << this->SliceRangeX[0] << " " << this->SliceRangeX[1] << "\n";
  os << indent << "SliceRangeY: "
     << this->SliceRangeY[0] << " " << this->SliceRangeY[1] << "\n";
  os << indent << "SliceRangeZ: "
     << this->SliceRangeZ[0] << " " << this->SliceRangeZ[1] << "\n";
  os << indent << "Stencil: " << this->GetStencil() << "\n";
  os << indent << "ActiveComponent: " << this->ActiveComponent << "\n";
}
