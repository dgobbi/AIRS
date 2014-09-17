/*=========================================================================

  Program:   Atamai Image Registration and Segmentation
  Module:    vtkImageConnectivityFilter.cxx

  Copyright (c) 2014 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkImageConnectivityFilter.h"

#include "vtkMath.h"
#include "vtkImageData.h"
#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkImageStencilData.h"
#include "vtkImageRegionIterator.h"
#include "vtkImageIterator.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkIdTypeArray.h"
#include "vtkIntArray.h"
#include "vtkImageStencilData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkTemplateAliasMacro.h"
#include "vtkTypeTraits.h"
#include "vtkSmartPointer.h"
#include "vtkVersion.h"

#include <vector>
#include <stack>
#include <algorithm>

vtkStandardNewMacro(vtkImageConnectivityFilter);

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageConnectivityFilter::vtkImageConnectivityFilter()
{
  this->LabelMode = SeedScalar;
  this->ExtractionMode = SeededRegions;

  this->ScalarRange[0] = 0.5;
  this->ScalarRange[1] = VTK_DOUBLE_MAX;

  this->SizeRange[0] = 1;
#if VTK_MAJOR_VERSION >= 6
  this->SizeRange[1] = VTK_ID_MAX;
#else
  this->SizeRange[1] = VTK_LARGE_ID;
#endif

  this->LabelConstantValue = 255;

  this->ActiveComponent = 0;

  this->LabelScalarType = VTK_UNSIGNED_CHAR;

  this->ExtractedRegionLabels = vtkIntArray::New();
  this->ExtractedRegionSizes = vtkIdTypeArray::New();
  this->ExtractedRegionSeedIds = vtkIdTypeArray::New();

  this->SetNumberOfInputPorts(3);
}

//----------------------------------------------------------------------------
vtkImageConnectivityFilter::~vtkImageConnectivityFilter()
{
  if (this->ExtractedRegionSizes)
    {
    this->ExtractedRegionSizes->Delete();
    }
  if (this->ExtractedRegionLabels)
    {
    this->ExtractedRegionLabels->Delete();
    }
  if (this->ExtractedRegionSeedIds)
    {
    this->ExtractedRegionSeedIds->Delete();
    }
}

//----------------------------------------------------------------------------
int vtkImageConnectivityFilter::FillInputPortInformation(
  int port, vtkInformation* info)
{
  if (port == 2)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }
  else if (port == 1)
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
void vtkImageConnectivityFilter::SetStencilConnection(
  vtkAlgorithmOutput *stencil)
{
  this->SetInputConnection(1, stencil);
}

//----------------------------------------------------------------------------
vtkAlgorithmOutput *vtkImageConnectivityFilter::GetStencilConnection()
{
  return this->GetInputConnection(1, 0);
}

//----------------------------------------------------------------------------
void vtkImageConnectivityFilter::SetStencilData(vtkImageStencilData *stencil)
{
#if VTK_MAJOR_VERSION >= 6
  this->SetInputData(1, stencil);
#else
  this->SetInput(1, stencil);
#endif
}

//----------------------------------------------------------------------------
void vtkImageConnectivityFilter::SetSeedConnection(
  vtkAlgorithmOutput *seeds)
{
  this->SetInputConnection(2, seeds);
}

//----------------------------------------------------------------------------
vtkAlgorithmOutput *vtkImageConnectivityFilter::GetSeedConnection()
{
  return this->GetInputConnection(2, 0);
}

//----------------------------------------------------------------------------
void vtkImageConnectivityFilter::SetSeedData(vtkPointSet *seeds)
{
#if VTK_MAJOR_VERSION >= 6
  this->SetInputData(2, seeds);
#else
  this->SetInput(2, seeds);
#endif
}

//----------------------------------------------------------------------------
const char *vtkImageConnectivityFilter::GetLabelScalarTypeAsString()
{
  const char *result = "Unknown";
  switch (this->LabelScalarType)
    {
    case VTK_UNSIGNED_CHAR:
      result = "UnsignedChar";
      break;
    case VTK_SHORT:
      result = "Short";
      break;
    case VTK_UNSIGNED_SHORT:
      result = "UnsignedShort";
      break;
    case VTK_INT:
      result = "Int";
      break;
    }
  return result;
}

//----------------------------------------------------------------------------
const char *vtkImageConnectivityFilter::GetLabelModeAsString()
{
  const char *result = "Unknown";
  switch (this->LabelMode)
    {
    case SeedScalar:
      result = "SeedScalar";
      break;
    case ConstantValue:
      result = "ConstantValue";
      break;
    case SizeRank:
      result = "SizeRank";
      break;
    }
  return result;
}

//----------------------------------------------------------------------------
const char *vtkImageConnectivityFilter::GetExtractionModeAsString()
{
  const char *result = "Unknown";
  switch (this->ExtractionMode)
    {
    case SeededRegions:
      result = "SeededRegions";
      break;
    case AllRegions:
      result = "AllRegions";
      break;
    case LargestRegion:
      result = "LargestRegion";
      break;
    }
  return result;
}

//----------------------------------------------------------------------------
vtkIdType vtkImageConnectivityFilter::GetNumberOfExtractedRegions()
{
  return this->ExtractedRegionLabels->GetNumberOfTuples();
}

//----------------------------------------------------------------------------
namespace {

// A private class for the connectivity algorithm
class vtkICF
{
public:
  // Simple struct that holds information about a region.
  struct Region;

  // A class that is a vector of regions.
  class RegionVector;

protected:
  // Simple class that holds a seed location and a scalar value.
  class Seed;

  // A functor to assist in comparing region sizes.
  struct CompareSize;

  // Create a bit mask from the provided stencil.
  static void MakeMaskFromStencil(
    unsigned char *maskPtr, vtkImageStencilData *stencil, int extent[6]);

  // Remove all but the largest region from the output image.
  template<class OT>
  static void PruneAllButLargest(
    vtkImageData *outData, OT *outPtr, vtkImageStencilData *stencil,
    int extent[6], const OT& value, vtkICF::RegionVector& regionInfo);

  // Remove the smallest region from the output image.
  // This is called when there are no labels left, i.e. when the label
  // value reaches the maximum allowed by the output data type.
  template<class OT>
  static void PruneSmallestRegion(
    vtkImageData *outData, OT *outPtr, vtkImageStencilData *stencil,
    int extent[6], vtkICF::RegionVector& regionInfo);

  // Remove all islands that aren't in the given range of sizes
  template<class OT>
  static void PruneBySize(
    vtkImageData *outData, OT *outPtr, vtkImageStencilData *stencil,
    int extent[6], vtkIdType sizeRange[2], vtkICF::RegionVector& regionInfo);

  // This is the function that grows a region from a seed.
  template<class IT, class OT>
  static vtkIdType Fill(
    IT *inPtr, vtkIdType inInc[3], IT low, IT high,
    OT *outPtr, vtkIdType outInc[3], int outLimits[6],
    unsigned char *maskPtr, int maxIdx[3],
    std::stack<vtkICF::Seed> &seedStack);

  // Add a region to the list of regions.
  template<class OT>
  static void AddRegion(
    vtkImageData *outData, OT *outPtr, vtkImageStencilData *stencil,
    int extent[6], vtkIdType sizeRange[2], vtkICF::RegionVector& regionInfo,
    vtkIdType voxelCount, vtkIdType regionId, int extractionMode);

  // Fill the ExtractedRegionSizes and ExtractedRegionLabels arrays.
  static void GenerateRegionArrays(
    vtkImageConnectivityFilter *self, vtkICF::RegionVector& regionInfo,
    vtkDataArray *seedScalars, int minLabel, int maxLabel);

  // Relabel the image, usually the last method called.
  template<class OT>
  static void Relabel(
    vtkImageData *outData, OT *outPtr,
    vtkImageStencilData *stencil, int extent[6],
    vtkIntArray *labelMap);

  // Sort the ExtractedRegionLabels array and the other arrays.
  static void SortRegionArrays(vtkImageConnectivityFilter *self);

  // Finalize the output
  template<class OT>
  static void Finish(
    vtkImageConnectivityFilter *self, vtkImageData *outData,
    OT *outPtr, vtkImageStencilData *stencil, int extent[6],
    vtkDataArray *seedScalars, vtkICF::RegionVector& regionInfo);

  // Execute method for when point seeds are provided.
  template <class IT, class OT>
  static void SeededExecute(
    vtkImageConnectivityFilter *self,
    vtkImageData *inData, vtkImageData *outData,
    vtkPointSet *seedData, vtkImageStencilData *stencil,
    IT *inPtr, OT *outPtr, unsigned char *maskPtr, int extent[6],
    IT srange[2], vtkICF::RegionVector& regionInfo, int id);

  // Execute method for when no seeds are provided.
  template <class IT, class OT>
  static void SeedlessExecute(
    vtkImageConnectivityFilter *self,
    vtkImageData *inData, vtkImageData *outData, vtkImageStencilData *stencil,
    IT *inPtr, OT *outPtr, unsigned char *maskPtr, int extent[6],
    IT srange[2], vtkICF::RegionVector& regionInfo, int id);

  // The main execute method.
  template <class IT, class OT>
  static void ExecuteInputOutput(
    vtkImageConnectivityFilter *self,
    vtkImageData *inData, vtkImageData *outData,
    vtkPointSet *seedData, vtkImageStencilData *stencil,
    IT *inPtr, IT srange[2], OT *outPtr, unsigned char *maskPtr,
    int extent[6], int id);

public:
  // Execute method templated only over input type (calls the above method).
  template <class IT>
  static void Execute(
    vtkImageConnectivityFilter *self,
    vtkImageData *inData, vtkImageData *outData,
    vtkPointSet *seedData, vtkImageStencilData *stencil,
    IT *inPtr, unsigned char *maskPtr, int extent[6], int id);

  // Utility method to find the intersection of two extents.
  // Returns false if the extents do not intersect.
  static bool IntersectExtents(
    const int extent1[6], const int extent2[6], int output[6]);
};

//----------------------------------------------------------------------------
// seed struct: structured coordinates and a value,
// the coords can be accessed with [] and the value with *
class vtkICF::Seed
{
public:
  Seed() {
    pos[0] = 0; pos[1] = 0; pos[2] = 0; value = 0; }

  Seed(int i, int j, int k) {
    pos[0] = i; pos[1] = j; pos[2] = k; value = 0; }

  Seed(int i, int j, int k, int v) {
    pos[0] = i; pos[1] = j; pos[2] = k; value = v; }

  Seed(const vtkICF::Seed &seed) {
    pos[0] = seed.pos[0];
    pos[1] = seed.pos[1];
    pos[2] = seed.pos[2];
    value = seed.value; };

  const int &operator[](int i) const { return pos[i]; }
  int &operator[](int i) { return pos[i]; }

  const int &operator*() const { return value; }
  int &operator*() { return value; }

  const vtkICF::Seed &operator = (const vtkICF::Seed seed) {
    pos[0] = seed.pos[0];
    pos[1] = seed.pos[1];
    pos[2] = seed.pos[2];
    value = seed.value;
    return *this; }

private:
  int pos[3];
  int value;
};

//----------------------------------------------------------------------------
// region struct: size and id
struct vtkICF::Region
{
  Region(vtkIdType s, vtkIdType i) : size(s), id(i) {}
  Region() : size(0), id(0) {}

  vtkIdType size;
  vtkIdType id;
};

//----------------------------------------------------------------------------
class vtkICF::RegionVector : public std::vector<vtkICF::Region>
{
public:
  typedef std::vector<vtkICF::Region>::iterator iterator;

  iterator smallest()
  {
    iterator small = begin() + 1;
    if (small != end())
      {
      for (iterator i = small + 1; i != end(); ++i)
        {
        if (i->size <= small->size)
          {
          small = i;
          }
        }
      }
    return small;
  }

  iterator largest()
  {
    iterator large = begin() + 1;
    if (large != end())
      {
      for (iterator i = large + 1; i != end(); ++i)
        {
        if (i->size > large->size)
          {
          large = i;
          }
        }
      }
    return large;
  }

};

//----------------------------------------------------------------------------
bool vtkICF::IntersectExtents(
  const int extent1[6], const int extent2[6], int output[6])
{
  bool rval = true;

  int i = 0;
  while (i < 6)
    {
    output[i] = (extent1[i] > extent2[i] ? extent1[i] : extent2[i]);
    i++;
    output[i] = (extent1[i] < extent2[i] ? extent1[i] : extent2[i]);
    rval &= (output[i-1] <= output[i]);
    i++;
    }

  return rval;
}

//----------------------------------------------------------------------------
void vtkICF::MakeMaskFromStencil(
  unsigned char *maskPtr, vtkImageStencilData *stencil, int extent[6])
{
  // The vtkImageRegionIterator requires a vtkImageData, so make
  // a fake vtkImageData for the stencil, it will contain only
  // one real voxel but its extent will be our desired extent
  vtkSmartPointer<vtkImageData> data =
    vtkSmartPointer<vtkImageData>::New();
  data->SetExtent(0,0,0,0,0,0);
#if VTK_MAJOR_VERSION >= 6
  data->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
#else
  data->SetScalarType(VTK_UNSIGNED_CHAR);
  data->AllocateScalars();
#endif
  data->SetExtent(extent);

  // for pointer math, get pointer to first voxel
  unsigned char *beginptr = static_cast<unsigned char *>(
    data->GetScalarPointer());

  // then set all bits that are outside the stencil region
  vtkImageRegionIterator<unsigned char> iter(data, stencil, extent);
  for (; !iter.IsAtEnd(); iter.NextSpan())
    {
    unsigned char *ptr = iter.BeginSpan();
    unsigned char *endptr = iter.EndSpan();
    if (!iter.IsInStencil())
      {
      for (; ptr != endptr; ptr++)
        {
        size_t offset = ptr - beginptr;
        unsigned char bitpos = static_cast<unsigned char>(offset) & 0x7;
        unsigned char *maskPtr1 = maskPtr + (offset >> 3);
        *maskPtr1 = *maskPtr1 ^ (1 << bitpos);
        }
      }
    }
}

//----------------------------------------------------------------------------
template<class OT>
void vtkICF::PruneAllButLargest(
  vtkImageData *outData, OT *outPtr, vtkImageStencilData *stencil,
  int extent[6], const OT& value, vtkICF::RegionVector& regionInfo)
{
  // clip the extent with the output extent
  int outExt[6];
  outData->GetExtent(outExt);
  if (!vtkICF::IntersectExtents(outExt, extent, outExt))
    {
    return;
    }

  // find the largest region
  vtkICF::RegionVector::iterator largest = regionInfo.largest();
  if (largest != regionInfo.end())
    {
    // get its position, remove all other regions from the list
    OT t = std::distance(regionInfo.begin(), largest);
    regionInfo[1] = *largest;
    regionInfo.erase(regionInfo.begin()+2, regionInfo.end());

    // remove all other regions from the output
    vtkImageRegionIterator<OT> iter(outData, stencil, outExt);
    for (; !iter.IsAtEnd(); iter.NextSpan())
      {
      if (iter.IsInStencil())
        {
        outPtr = iter.BeginSpan();
        OT *endPtr = iter.EndSpan();
        for (; outPtr != endPtr; ++outPtr)
          {
          OT v = *outPtr;
          if (v == t)
            {
            *outPtr = value;
            }
          else if (v != 0)
            {
            *outPtr = 0;
            }
          }
        }
      }
    }
}

//----------------------------------------------------------------------------
template<class OT>
void vtkICF::PruneSmallestRegion(
  vtkImageData *outData, OT *outPtr, vtkImageStencilData *stencil,
  int extent[6], vtkICF::RegionVector& regionInfo)
{
  // clip the extent with the output extent
  int outExt[6];
  outData->GetExtent(outExt);
  if (!vtkICF::IntersectExtents(outExt, extent, outExt))
    {
    return;
    }

  // find the smallest region
  vtkICF::RegionVector::iterator smallest = regionInfo.smallest();
  if (smallest != regionInfo.end())
    {
    // get the index to the smallest value and remove it
    OT t = std::distance(regionInfo.begin(), smallest);
    regionInfo.erase(smallest);

    // remove the corresponding region from the output
    vtkImageRegionIterator<OT> iter(outData, stencil, outExt);
    for (; !iter.IsAtEnd(); iter.NextSpan())
      {
      if (iter.IsInStencil())
        {
        outPtr = iter.BeginSpan();
        OT *endPtr = iter.EndSpan();
        for (; outPtr != endPtr; ++outPtr)
          {
          OT v = *outPtr;
          if (v == t)
            {
            *outPtr = 0;
            }
          else if (v > t)
            {
            *outPtr = v-1;
            }
          }
        }
      }
    }
}

//----------------------------------------------------------------------------
template<class OT>
void vtkICF::PruneBySize(
  vtkImageData *outData, OT *outPtr, vtkImageStencilData *stencil,
  int extent[6], vtkIdType sizeRange[2], vtkICF::RegionVector& regionInfo)
{
  // find all the regions in the allowed size range
  size_t n = regionInfo.size();
  std::vector<OT> newlabels(n);
  newlabels[0] = 0;
  size_t m = 1;
  for (size_t i = 1; i < n; i++)
    {
    size_t l = 0;
    vtkIdType s = regionInfo[i].size;
    if (s >= sizeRange[0] && s <= sizeRange[1])
      {
      l = m++;
      if (i != l)
        {
        regionInfo[l] = regionInfo[i];
        }
      }
    newlabels[i] = static_cast<OT>(l);
    }

  // were any regions outside of the range?
  if (m < n)
    {
    // resize regionInfo
    regionInfo.resize(m);

    // clip the extent with the output extent
    int outExt[6];
    outData->GetExtent(outExt);
    if (!vtkICF::IntersectExtents(outExt, extent, outExt))
      {
      return;
      }

    // remove the corresponding regions from the output
    vtkImageRegionIterator<OT> iter(outData, stencil, outExt);
    for (; !iter.IsAtEnd(); iter.NextSpan())
      {
      if (iter.IsInStencil())
        {
        outPtr = iter.BeginSpan();
        OT *endPtr = iter.EndSpan();
        for (; outPtr != endPtr; ++outPtr)
          {
          OT v = *outPtr;
          if (v != 0)
            {
            *outPtr = newlabels[v];
            }
          }
        }
      }
    }
}

//----------------------------------------------------------------------------
// Perform a flood fill for each given seed.
template<class IT, class OT>
vtkIdType vtkICF::Fill(
  IT *inPtr, vtkIdType inInc[3], IT low, IT high,
  OT *outPtr, vtkIdType outInc[3], int outLimits[6],
  unsigned char *maskPtr, int maxIdx[3],
  std::stack<vtkICF::Seed> &seedStack)
{
  vtkIdType counter = 0;

  while (!seedStack.empty())
    {
    // get the seed at the top of the stack
    vtkICF::Seed seed = seedStack.top();
    seedStack.pop();

    // get the value of the input voxel at this position
    vtkIdType inOffset = (seed[0]*inInc[0] +
                          seed[1]*inInc[1] +
                          seed[2]*inInc[2]);
    IT value = inPtr[inOffset];

    // make sure the input is within the threshold
    if (value < low || value > high)
      {
      continue;
      }

    // get the pointer into the bitmask
    unsigned char bit = 1 << static_cast<unsigned char>(inOffset & 0x7);
    unsigned char *maskPtr1 = maskPtr + (inOffset >> 3);

    // if already colored, skip
    if ((*maskPtr1 & bit) != 0)
      {
      continue;
      }

    // paint the mask and count the voxel
    *maskPtr1 ^= bit;
    counter++;

    // set the output if within output extent
    if (seed[0] >= outLimits[0] && seed[0] <= outLimits[1] &&
        seed[1] >= outLimits[2] && seed[1] <= outLimits[3] &&
        seed[2] >= outLimits[4] && seed[2] <= outLimits[5])
      {
      // get the pointer to the seed position in the mask
      vtkIdType outOffset = ((seed[0] - outLimits[0])*outInc[0] +
                             (seed[1] - outLimits[2])*outInc[1] +
                             (seed[2] - outLimits[4])*outInc[2]);

      outPtr[outOffset] = static_cast<OT>(*seed);
      }

    // push the new seeds for the six neighbors, make sure offsets in X are
    // pushed last so that they will be popped first (we want to raster X,
    // Y, and Z in that order).
    int i = 3;
    do
      {
      i--;
      if (seed[i] > 0)
        {
        seed[i]--;
        seedStack.push(seed);
        seed[i]++;
        }
      if (seed[i] < maxIdx[i])
        {
        seed[i]++;
        seedStack.push(seed);
        seed[i]--;
        }
      }
    while (i != 0);
    }

  return counter;
}

//----------------------------------------------------------------------------
template<class OT>
void vtkICF::AddRegion(
  vtkImageData *outData, OT *outPtr, vtkImageStencilData *stencil,
  int extent[6], vtkIdType sizeRange[2], vtkICF::RegionVector& regionInfo,
  vtkIdType voxelCount, vtkIdType regionId, int extractionMode)
{
  regionInfo.push_back(vtkICF::Region(voxelCount, regionId));
  // check if the label value has reached its maximum, and if so,
  // remove some of the regions
  if (regionInfo.size() > static_cast<size_t>(vtkTypeTraits<OT>::Max()))
    {
    vtkICF::PruneBySize(
      outData, outPtr, stencil, extent, sizeRange, regionInfo);

    // if that didn't remove anything, try these:
    if (regionInfo.size() > static_cast<size_t>(vtkTypeTraits<OT>::Max()))
      {
      if (extractionMode == vtkImageConnectivityFilter::LargestRegion)
        {
        OT label = 1;
        vtkICF::PruneAllButLargest(
          outData, outPtr, stencil, extent, label, regionInfo);
        }
      else
        {
        vtkICF::PruneSmallestRegion(
          outData, outPtr, stencil, extent, regionInfo);
        }
      }
    }
}

//----------------------------------------------------------------------------
// a functor to sort region indices by region size
struct vtkICF::CompareSize
{
  CompareSize(vtkICF::RegionVector &r)
    : Regions(&r) {}

  vtkICF::RegionVector *Regions;

  bool operator()(vtkIdType x, vtkIdType y)
  {
    return ((*Regions)[x].size > (*Regions)[y].size);
  }
};

//----------------------------------------------------------------------------
void vtkICF::GenerateRegionArrays(
  vtkImageConnectivityFilter *self, vtkICF::RegionVector& regionInfo,
  vtkDataArray *seedScalars, int minLabel, int maxLabel)
{
  // clamp the default label value to the range of the output data type
  int constantLabel = self->GetLabelConstantValue();
  constantLabel = (constantLabel > minLabel ? constantLabel : minLabel);
  constantLabel = (constantLabel < maxLabel ? constantLabel : maxLabel);

  // build the arrays
  vtkIdTypeArray *sizes = self->GetExtractedRegionSizes();
  vtkIdTypeArray *ids = self->GetExtractedRegionSeedIds();
  vtkIntArray *labels = self->GetExtractedRegionLabels();

  if (regionInfo.size() == 1)
    {
    // only background is present, there are no connected regions
    sizes->Reset();
    ids->Reset();
    labels->Reset();
    }
  else if (self->GetExtractionMode() ==
           vtkImageConnectivityFilter::LargestRegion)
    {
    // only one region (the largest) will be output
    sizes->SetNumberOfValues(1);
    ids->SetNumberOfValues(1);
    labels->SetNumberOfValues(1);

    // get info for the largest region
    vtkICF::RegionVector::iterator largest = regionInfo.largest();

    // default label value is 1
    int label = 1;

    // check which label mode was selected
    switch (self->GetLabelMode())
      {
      // use the scalars of the seed points as labels
      case vtkImageConnectivityFilter::SeedScalar:
        if (seedScalars)
          {
          label = constantLabel;
          // get the label from the seed scalars
          if (largest->id >= 0)
            {
            double s = seedScalars->GetTuple1(largest->id);
            s = (s > minLabel ? s : minLabel);
            s = (s < maxLabel ? s : maxLabel);
            label = vtkMath::Floor(s + 0.5);
            }
          }
        break;
      // use the specified ConstantValue for all regions
      case vtkImageConnectivityFilter::ConstantValue:
        label = constantLabel;
        break;
      }

    // create the arrays for the single region present in the output
    sizes->SetValue(0, largest->size);
    ids->SetValue(0, largest->id);
    labels->SetValue(0, label);
    }
  else
    {
    // multiple regions might be present in the output (we subtract
    // one because the background doesn't count as a region)
    vtkIdType n = static_cast<vtkIdType>(regionInfo.size()) - 1;
    sizes->SetNumberOfValues(n);
    ids->SetNumberOfValues(n);
    labels->SetNumberOfValues(n);

    // build the arrays (this part is easy!)
    for (vtkIdType i = 0; i < n; i++)
      {
      sizes->SetValue(i, regionInfo[i+1].size);
      ids->SetValue(i, regionInfo[i+1].id);
      labels->SetValue(i, i+1);
      }

    // some label modes require additional actions to be done
    switch (self->GetLabelMode())
      {
      // change the labels to match the scalars for the seed points
      case vtkImageConnectivityFilter::SeedScalar:
        if (seedScalars)
          {
          for (vtkIdType i = 0; i < n; i++)
            {
            int label = constantLabel;
            vtkIdType j = regionInfo[i+1].id;
            if (j >= 0)
              {
              double s = seedScalars->GetTuple1(j);
              s = (s > minLabel ? s : minLabel);
              s = (s < maxLabel ? s : maxLabel);
              label = vtkMath::Floor(s + 0.5);
              }
            labels->SetValue(i, label);
            }
          }
        break;
      // order the labels by the size rank of the regions
      case vtkImageConnectivityFilter::SizeRank:
        {
        vtkICF::CompareSize cmpfunc(regionInfo);
        int *la = labels->GetPointer(0);
        std::stable_sort(la, la+n, cmpfunc);
        }
        break;
      // set all labels to the same value
      case vtkImageConnectivityFilter::ConstantValue:
        for (vtkIdType i = 0; i < n; i++)
          {
          labels->SetValue(i, constantLabel);
          }
        break;
      }
    }
}

//----------------------------------------------------------------------------
// generate the output image
template<class OT>
void vtkICF::Relabel(
  vtkImageData *outData, OT *outPtr,
  vtkImageStencilData *stencil, int extent[6],
  vtkIntArray *labelMap)
{
  // clip the extent with the output extent
  int outExt[6];
  outData->GetExtent(outExt);
  if (!vtkICF::IntersectExtents(outExt, extent, outExt))
    {
    return;
    }

  vtkImageRegionIterator<OT> iter(outData, stencil, outExt);

  // loop through the output voxels and change the "region id" value
  // stored in the voxel into a "region label" value.
  for (; !iter.IsAtEnd(); iter.NextSpan())
    {
    outPtr = iter.BeginSpan();
    OT *outEnd = iter.EndSpan();

    if (iter.IsInStencil())
      {
      for (; outPtr != outEnd; outPtr++)
        {
        OT v = *outPtr;
        if (v > 0)
          {
          *outPtr = labelMap->GetValue(v-1);
          }
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkICF::SortRegionArrays(
  vtkImageConnectivityFilter *self)
{
  vtkIdTypeArray *sizes = self->GetExtractedRegionSizes();
  vtkIdTypeArray *ids = self->GetExtractedRegionSeedIds();
  vtkIntArray *labels = self->GetExtractedRegionLabels();

  vtkIdType *sizePtr = sizes->GetPointer(0);
  vtkIdType *idPtr = ids->GetPointer(0);
  int *labelPtr = labels->GetPointer(0);

  vtkIdType n = labels->GetNumberOfTuples();

  // only SizeRank needs re-sorting
  if (self->GetLabelMode() == vtkImageConnectivityFilter::SizeRank)
    {
    std::vector<vtkIdType> sizeVector(sizePtr, sizePtr + n);
    std::vector<vtkIdType> idVector(idPtr, idPtr + n);
    for (vtkIdType i = 0; i < n; i++)
      {
      int j = labelPtr[i] - 1;
      sizePtr[i] = sizeVector[j];
      idPtr[i] = idVector[j];
      labelPtr[i] = static_cast<int>(i + 1);
      }
    }
}

//----------------------------------------------------------------------------
template<class OT>
void vtkICF::Finish(
  vtkImageConnectivityFilter *self, vtkImageData *outData,
  OT *outPtr, vtkImageStencilData *stencil, int extent[6],
  vtkDataArray *seedScalars, vtkICF::RegionVector& regionInfo)
{
  // Get the execution parameters
  int labelMode = self->GetLabelMode();
  int extractionMode = self->GetExtractionMode();
  vtkIdType sizeRange[2];
  self->GetSizeRange(sizeRange);

  // get only the regions in the requested range of sizes
  vtkICF::PruneBySize(
    outData, outPtr, stencil, extent, sizeRange, regionInfo);

  // create the three region info arrays
  vtkICF::GenerateRegionArrays(
    self, regionInfo, seedScalars,
    vtkTypeTraits<OT>::Min(), vtkTypeTraits<OT>::Max());

  vtkIntArray *labelArray = self->GetExtractedRegionLabels();
  if (labelArray->GetNumberOfTuples() > 0)
    {
    // do the extraction and final labeling
    if (extractionMode == vtkImageConnectivityFilter::LargestRegion)
      {
      OT label = static_cast<OT>(labelArray->GetValue(0));
      vtkICF::PruneAllButLargest(
        outData, outPtr, stencil, extent, label, regionInfo);
      }
    else if (labelMode != vtkImageConnectivityFilter::SeedScalar ||
             seedScalars != 0)
      {
      // this is done unless labelMode == SeedScalar and seedScalars == 0
      vtkICF::Relabel(
        outData, outPtr, stencil, extent, labelArray);
      }

    // sort the three region info arrays (must be done after Relabel)
    vtkICF::SortRegionArrays(self);
    }

}

//----------------------------------------------------------------------------
template <class IT, class OT>
void vtkICF::SeededExecute(
  vtkImageConnectivityFilter *self,
  vtkImageData *inData, vtkImageData *outData,
  vtkPointSet *seedData, vtkImageStencilData *stencil,
  IT *inPtr, OT *outPtr, unsigned char *maskPtr, int extent[6],
  IT srange[2], vtkICF::RegionVector& regionInfo, int)
{
  // Get execution parameters
  int extractionMode = self->GetExtractionMode();
  vtkIdType sizeRange[2];
  self->GetSizeRange(sizeRange);

  vtkIdType inInc[3];
  vtkIdType outInc[3];
  inData->GetIncrements(inInc);
  outData->GetIncrements(outInc);

  double spacing[3];
  double origin[3];
  inData->GetOrigin(origin);
  inData->GetSpacing(spacing);

  // Indexing goes from 0 to maxIdX
  int maxIdx[3];
  maxIdx[0] = extent[1] - extent[0];
  maxIdx[1] = extent[3] - extent[2];
  maxIdx[2] = extent[5] - extent[4];

  // Get the limits for the output data
  int outLimits[6];
  outData->GetExtent(outLimits);
  outLimits[0] -= extent[0];
  outLimits[1] -= extent[0];
  outLimits[2] -= extent[2];
  outLimits[3] -= extent[2];
  outLimits[4] -= extent[4];
  outLimits[5] -= extent[4];

  // label consecutively, starting at 1
  OT label = 1;

  std::stack<vtkICF::Seed> seedStack;

  vtkPoints *points = seedData->GetPoints();
  vtkIdType nPoints = points->GetNumberOfPoints();

  for (vtkIdType i = 0; i < nPoints; i++)
    {
    double point[3];
    points->GetPoint(i, point);
    int idx[3];
    bool outOfBounds = false;

    // convert point from data coords to image index
    for (int j = 0; j < 3; j++)
      {
      idx[j] = vtkMath::Floor((point[j] - origin[j])/spacing[j] + 0.5);
      idx[j] -= extent[2*j];
      outOfBounds |= (idx[j] < 0 || idx[j] > maxIdx[j]);
      }

    if (outOfBounds)
      {
      continue;
      }

    seedStack.push(vtkICF::Seed(idx[0], idx[1], idx[2], label));

    // find all voxels that are connected to the seed
    vtkIdType voxelCount = vtkICF::Fill(
      inPtr, inInc, srange[0], srange[1],
      outPtr, outInc, outLimits, maskPtr, maxIdx,
      seedStack);

    if (voxelCount != 0)
      {
      vtkICF::AddRegion(
        outData, outPtr, stencil, extent, sizeRange, regionInfo,
        voxelCount, i, extractionMode);
      label = static_cast<OT>(regionInfo.size());
      }
    }
}

//----------------------------------------------------------------------------
template <class IT, class OT>
void vtkICF::SeedlessExecute(
  vtkImageConnectivityFilter *self,
  vtkImageData *inData, vtkImageData *outData, vtkImageStencilData *stencil,
  IT *inPtr, OT *outPtr, unsigned char *maskPtr, int extent[6],
  IT srange[2], vtkICF::RegionVector& regionInfo, int id)
{
  // Get execution parameters
  int extractionMode = self->GetExtractionMode();
  vtkIdType sizeRange[2];
  self->GetSizeRange(sizeRange);

  vtkIdType inInc[3];
  vtkIdType outInc[3];
  inData->GetIncrements(inInc);
  outData->GetIncrements(outInc);

  // Indexing goes from 0 to maxIdX
  int maxIdx[3];
  maxIdx[0] = extent[1] - extent[0];
  maxIdx[1] = extent[3] - extent[2];
  maxIdx[2] = extent[5] - extent[4];

  // Get the limits for the output data
  int outLimits[6];
  outData->GetExtent(outLimits);
  outLimits[0] -= extent[0];
  outLimits[1] -= extent[0];
  outLimits[2] -= extent[2];
  outLimits[3] -= extent[2];
  outLimits[4] -= extent[4];
  outLimits[5] -= extent[4];

  std::stack<vtkICF::Seed> seedStack;

  vtkImageRegionIterator<IT> iter(inData, stencil, extent, self, id);
  for (; !iter.IsAtEnd(); iter.NextSpan())
    {
    if (iter.IsInStencil())
      {
      int xIdx = iter.GetIndexX() - extent[0];
      int yIdx = iter.GetIndexY() - extent[2];
      int zIdx = iter.GetIndexZ() - extent[4];

      IT *inPtr0 = iter.BeginSpan();
      IT *inPtrEnd = iter.EndSpan();

      for (; inPtr0 != inPtrEnd; xIdx++, inPtr0 += inInc[0])
        {
        if (*inPtr0 < srange[0] || *inPtr0 > srange[1])
          {
          continue;
          }

        OT label = static_cast<OT>(regionInfo.size());
        seedStack.push(vtkICF::Seed(xIdx, yIdx, zIdx, label));

        // find all voxels that are connected to the seed
        vtkIdType voxelCount = vtkICF::Fill(
          inPtr, inInc, srange[0], srange[1],
          outPtr, outInc, outLimits, maskPtr, maxIdx,
          seedStack);

        if (voxelCount != 0)
          {
          if (voxelCount == 1 &&
              static_cast<OT>(regionInfo.size()) == vtkTypeTraits<OT>::Max())
            {
            // smallest region is definitely the one we just added
            vtkIdType outOffset = (xIdx*outInc[0] +
                                   yIdx*outInc[1] +
                                   zIdx*outInc[2]);
            outPtr[outOffset] = 0;
            }
          else
            {
            vtkICF::AddRegion(
              outData, outPtr, stencil, extent, sizeRange, regionInfo,
              voxelCount, -1, extractionMode);
            }
          }
        }
      }
    }
}

//----------------------------------------------------------------------------
// This templated function executes the filter for any type of data.
template <class IT, class OT>
void vtkICF::ExecuteInputOutput(
  vtkImageConnectivityFilter *self,
  vtkImageData *inData, vtkImageData *outData, vtkPointSet *seedData,
  vtkImageStencilData *stencil, IT *inPtr, IT srange[2],
  OT *outPtr, unsigned char *maskPtr, int extent[6], int id)
{
  // push the "background" onto the region vector
  vtkICF::RegionVector regionInfo;
  regionInfo.push_back(vtkICF::Region(0, 0));

  // execution depends on how regions are seeded
  vtkDataArray *seedScalars = 0;
  if (seedData)
    {
    seedScalars = seedData->GetPointData()->GetScalars();
    vtkICF::SeededExecute(
      self, inData, outData, seedData, stencil, inPtr, outPtr, maskPtr,
      extent, srange, regionInfo, id);
    }

  // if no seeds, or if AllRegions selected, search for all regions
  int extractionMode = self->GetExtractionMode();
  if (!seedData ||
      extractionMode == vtkImageConnectivityFilter::AllRegions)
    {
    vtkICF::SeedlessExecute(
      self, inData, outData, stencil, inPtr, outPtr, maskPtr, extent,
      srange, regionInfo, id);
    }

  // do final relabelling and other bookkeeping
  vtkICF::Finish(
    self, outData, outPtr, stencil, extent, seedScalars, regionInfo);
}

//----------------------------------------------------------------------------
template <class IT>
void vtkICF::Execute(
  vtkImageConnectivityFilter *self,
  vtkImageData *inData, vtkImageData *outData,
  vtkPointSet *seedData, vtkImageStencilData *stencil,
  IT *inPtr, unsigned char *maskPtr, int extent[6], int id)
{
  // Get active component (only one component is thresholded)
  int nComponents = inData->GetNumberOfScalarComponents();
  int activeComponent = self->GetActiveComponent();
  if (activeComponent < 0 || activeComponent > nComponents)
    {
    activeComponent = 0;
    }

  // Get the scalar range clamped to the input type range
  double drange[2];
  self->GetScalarRange(drange);
  IT srange[2];
  srange[0] = vtkTypeTraits<IT>::Min();
  srange[1] = vtkTypeTraits<IT>::Max();
  if (drange[0] > static_cast<double>(srange[1]))
    {
    srange[0] = srange[1];
    }
  else if (drange[0] > static_cast<double>(srange[0]))
    {
    srange[0] = static_cast<IT>(drange[0]);
    }
  if (drange[1] < static_cast<double>(srange[0]))
    {
    srange[1] = srange[0];
    }
  else if (drange[1] < static_cast<double>(srange[1]))
    {
    srange[1] = static_cast<IT>(drange[1]);
    }

  // Get pointer for input extent
  inPtr = static_cast<IT *>(inData->GetScalarPointerForExtent(extent));
  // Get pointer for entire output extent
  void *outPtr = outData->GetScalarPointer();

  // Adjust pointer to use the active component
  inPtr += activeComponent;

  // pre-set all mask voxels that are outside the stencil
  vtkICF::MakeMaskFromStencil(maskPtr, stencil, extent);

  switch (outData->GetScalarType())
    {
    case VTK_UNSIGNED_CHAR:
      vtkICF::ExecuteInputOutput(
        self, inData, outData, seedData, stencil, inPtr, srange,
        static_cast<unsigned char *>(outPtr), maskPtr, extent, id);
      break;

    case VTK_SHORT:
      vtkICF::ExecuteInputOutput(
        self, inData, outData, seedData, stencil, inPtr, srange,
        static_cast<short *>(outPtr), maskPtr, extent, id);
      break;

    case VTK_UNSIGNED_SHORT:
      vtkICF::ExecuteInputOutput(
        self, inData, outData, seedData, stencil, inPtr, srange,
        static_cast<unsigned short *>(outPtr), maskPtr, extent, id);
      break;

    case VTK_INT:
      vtkICF::ExecuteInputOutput(
        self, inData, outData, seedData, stencil, inPtr, srange,
        static_cast<int *>(outPtr), maskPtr, extent, id);
      break;
    }
}

} // end anonymous namespace

//----------------------------------------------------------------------------
int vtkImageConnectivityFilter::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkDataObject::SetPointDataActiveScalarInfo(
    outInfo, this->LabelScalarType, 1);

  return 1;
}

//----------------------------------------------------------------------------
int vtkImageConnectivityFilter::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector))
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *stencilInfo = inputVector[1]->GetInformationObject(0);

  int extent[6];
  inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent, 6);
  if (stencilInfo)
    {
    stencilInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
                     extent, 6);
    }

  return 1;
}

//----------------------------------------------------------------------------
int vtkImageConnectivityFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *stencilInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation *seedInfo = inputVector[2]->GetInformationObject(0);

  vtkImageData* outData = static_cast<vtkImageData *>(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkImageData* inData = static_cast<vtkImageData *>(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkPointSet* seedData = 0;
  if (seedInfo)
    {
    seedData = static_cast<vtkPointSet *>(
      seedInfo->Get(vtkDataObject::DATA_OBJECT()));
    }

  vtkImageStencilData* stencil = 0;
  if (stencilInfo)
    {
    stencil = static_cast<vtkImageStencilData *>(
      stencilInfo->Get(vtkDataObject::DATA_OBJECT()));
    }

  int outExt[6];
  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), outExt);
#if VTK_MAJOR_VERSION >= 6
  this->AllocateOutputData(outData, outInfo, outExt);
#else
  this->AllocateOutputData(outData, outExt);
#endif

  void *outPtr = outData->GetScalarPointerForExtent(outExt);

  // clear the output
  size_t size = (outExt[1] - outExt[0] + 1);
  size *= (outExt[3] - outExt[2] + 1);
  size *= (outExt[5] - outExt[4] + 1);
  size_t byteSize = size*outData->GetScalarSize();
  memset(outPtr, 0, byteSize);

  // we need all the voxels that might be connected to the seed
  int extent[6];
  inData->GetExtent(extent);

  // voxels outside the stencil extent can be excluded
  if (stencil)
    {
    int stencilExtent[6];
    stencil->GetExtent(stencilExtent);
    if (!vtkICF::IntersectExtents(extent, stencilExtent, extent))
      {
      // if stencil doesn't overlap the input, return
      return 1;
      }
    }

  int outScalarType = outData->GetScalarType();
  if (outScalarType != VTK_UNSIGNED_CHAR &&
      outScalarType != VTK_SHORT &&
      outScalarType != VTK_UNSIGNED_SHORT &&
      outScalarType != VTK_INT)
    {
    vtkErrorMacro("Execute: Output ScalarType is "
                  << outData->GetScalarType()
                  << ", but it must be one of VTK_UNSIGNED_CHAR, "
                  "VTK_SHORT, VTK_UNSIGNED_SHORT, or VTK_INT");
    return 0;
    }

  // create and clear the image bitmask (each bit is a voxel)
  size = (extent[1] - extent[0] + 1);
  size *= (extent[3] - extent[2] + 1);
  size *= (extent[5] - extent[4] + 1);
  byteSize = (size + 7) / 8;
  unsigned char *mask = new unsigned char [byteSize];
  memset(mask, 0, byteSize);

  // get scalar pointers
  void *inPtr = inData->GetScalarPointerForExtent(extent);

  int id = 0; // not multi-threaded, thread id is always zero

  int rval = 1;
  switch (inData->GetScalarType())
    {
    vtkTemplateAliasMacro(
      vtkICF::Execute(
        this, inData, outData, seedData, stencil,
        static_cast<VTK_TT *>(inPtr), mask, extent, id));

    default:
      vtkErrorMacro(<< "Execute: Unknown input ScalarType");
      rval = 0;
    }

  delete [] mask;

  return rval;
}

//----------------------------------------------------------------------------
void vtkImageConnectivityFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "LabelScalarType: "
     << this->GetLabelScalarTypeAsString() << "\n";

  os << indent << "LabelMode: "
     << this->GetLabelModeAsString() << "\n";

  os << indent << "ExtractionMode: "
     << this->GetExtractionModeAsString() << "\n";

  os << indent << "LabelConstantValue: "
     << this->LabelConstantValue << "\n";

  os << indent << "NumberOfExtractedRegions: "
     << this->GetNumberOfExtractedRegions() << "\n";

  os << indent << "ExtractedRegionLabels: "
     << this->ExtractedRegionLabels << "\n";

  os << indent << "ExtractedRegionSizes: "
     << this->ExtractedRegionSizes << "\n";

  os << indent << "ExtractedRegionSeedIds: "
     << this->ExtractedRegionSeedIds << "\n";

  os << indent << "ScalarRange: "
     << this->ScalarRange[0] << " " << this->ScalarRange[1] << "\n";

  os << indent << "SizeRange: "
     << this->SizeRange[0] << " " << this->SizeRange[1] << "\n";

  os << indent << "ActiveComponent: "
     << this->ActiveComponent << "\n";

  os << indent << "SeedConnection: "
     << this->GetSeedConnection() << "\n";

  os << indent << "StencilConnection: "
     << this->GetStencilConnection() << "\n";
}
