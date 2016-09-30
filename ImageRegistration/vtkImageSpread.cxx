/*=========================================================================

  Module: vtkImageSpread.cxx

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkImageSpread.h"

#include "vtkMath.h"
#include "vtkImageData.h"
#include "vtkImageStencilData.h"
#include "vtkImageStencilIterator.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkTemplateAliasMacro.h"
#include "vtkSmartPointer.h"
#include "vtkVersion.h"

vtkStandardNewMacro(vtkImageSpread);

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageSpread::vtkImageSpread()
{
  this->NumberOfIterations = 1;

  this->SetNumberOfInputPorts(2);
}

//----------------------------------------------------------------------------
vtkImageSpread::~vtkImageSpread()
{
}

//----------------------------------------------------------------------------
int vtkImageSpread::FillInputPortInformation(
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
void vtkImageSpread::SetStencilConnection(vtkAlgorithmOutput *stencil)
{
  this->SetInputConnection(1, stencil);
}

//----------------------------------------------------------------------------
vtkAlgorithmOutput *vtkImageSpread::GetStencilConnection()
{
  return this->GetInputConnection(1, 0);
}

//----------------------------------------------------------------------------
void vtkImageSpread::SetStencilData(vtkImageStencilData *stencil)
{
#if VTK_MAJOR_VERSION >= 6
  this->SetInputData(1, stencil);
#else
  this->SetInput(1, stencil);
#endif
}

//----------------------------------------------------------------------------
namespace {

// Internal methods for the algorithm
struct vtkISF
{
  // Main execution method
  template<class IT>
  static void Execute(
    vtkImageSpread *self, vtkImageData *inData, IT *inPtr,
    vtkImageData *outData, IT *outPtr, unsigned char *maskPtr,
    vtkImageStencilData *stencil, int extent[6], int iterations);

  // Utility method to find the intersection of two extents.
  // Returns false if the extents do not intersect.
  static bool IntersectExtents(
    const int extent1[6], const int extent2[6], int output[6]);

  // Rounding
  template<class T>
  static void Convert(double val, T *o)
  {
    *o = static_cast<T>(vtkMath::Round(val));
  }

  // No rounding for double
  static void Convert(double val, double *o)
  {
    *o = val;
  }

  // No rounding for float
  static void Convert(double val, float *o)
  {
    *o = val;
  }
};

//----------------------------------------------------------------------------
bool vtkISF::IntersectExtents(
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
template<class T>
void vtkISF::Execute(
  vtkImageSpread *, vtkImageData *inData, T *,
  vtkImageData *outData, T *outPtr, unsigned char *maskPtr,
  vtkImageStencilData *stencil, int extent[6], int iterations)
{
  int numComponents = outData->GetNumberOfScalarComponents();

  // copy everything inside the stencil to the output, create the mask
  unsigned char *maskPtrTmp = maskPtr;
  T *outPtrTmp = outPtr;
  vtkImageStencilIterator<T> iter(inData, stencil, extent);
  for (; !iter.IsAtEnd(); iter.NextSpan())
    {
    T *inPtr = iter.BeginSpan();
    T *inPtrEnd = iter.EndSpan();
    if (iter.IsInStencil())
      {
      while (inPtr != inPtrEnd)
        {
        *maskPtrTmp++ = 1;
        *outPtrTmp = *inPtr;
        for (int i = 1; i < numComponents; i++)
          {
          outPtrTmp[i] = inPtr[i];
          }
        outPtrTmp += numComponents;
        inPtr += numComponents;
        }
      }
    else
      {
      while (inPtr != inPtrEnd)
        {
        *maskPtrTmp++ = 0;
        *outPtrTmp = 0;
        for (int i = 1; i < numComponents; i++)
          {
          outPtrTmp[i] = 0;
          }
        outPtrTmp += numComponents;
        inPtr += numComponents;
        }
      }
    }

  // set "outside" voxels to the average of their neighbors
  int maxX = extent[1] - extent[0];
  int maxY = extent[3] - extent[2];
  int maxZ = extent[5] - extent[4];
  vtkIdType maskIncY = (maxX + 1);
  vtkIdType maskIncZ = maskIncY*(maxY + 1);
  vtkIdType outIncY = maskIncY*numComponents;
  vtkIdType outIncZ = maskIncZ*numComponents;

  double *pixel = new double[numComponents];

  unsigned char lastMask = 1;
  unsigned char newMask = 2;
  bool hit = false;

  for (int k = 0; k < iterations; k++)
    {
    outPtrTmp = outPtr;
    maskPtrTmp = maskPtr;

    for (int idZ = 0; idZ <= maxZ; idZ++)
      {
      for (int idY = 0; idY <= maxY; idY++)
        {
        for (int idX = 0; idX <= maxX; idX++)
          {
          int count = 0;
          *pixel = 0.0;
          for (int i = 1; i < numComponents; i++)
            {
            pixel[i] = 0.0;
            }

          if (*maskPtrTmp == 0)
            {
            if (idX > 0 && maskPtrTmp[-1] == lastMask)
              {
              count++;
              T *ptr = outPtrTmp - numComponents;
              *pixel += *ptr;
              for (int i = 1; i < numComponents; i++)
                {
                pixel[i] += ptr[i];
                }
              }
            if (idX < maxX && maskPtrTmp[1] == lastMask)
              {
              count++;
              T *ptr = outPtrTmp + numComponents;
              *pixel += *ptr;
              for (int i = 1; i < numComponents; i++)
                {
                pixel[i] += ptr[i];
                }
              }
            if (idY > 0 && maskPtrTmp[-maskIncY] == lastMask)
              {
              count++;
              T *ptr = outPtrTmp - outIncY;
              *pixel += *ptr;
              for (int i = 1; i < numComponents; i++)
                {
                pixel[i] += ptr[i];
                }
              }
            if (idY < maxY && maskPtrTmp[maskIncY] == lastMask)
              {
              count++;
              T *ptr = outPtrTmp + outIncY;
              *pixel += *ptr;
              for (int i = 1; i < numComponents; i++)
                {
                pixel[i] += ptr[i];
                }
              }
            if (idZ > 0 && maskPtrTmp[-maskIncZ] == lastMask)
              {
              count++;
              T *ptr = outPtrTmp - outIncZ;
              *pixel += *ptr;
              for (int i = 1; i < numComponents; i++)
                {
                pixel[i] += ptr[i];
                }
              }
            if (idZ < maxZ && maskPtrTmp[maskIncZ] == lastMask)
              {
              count++;
              T *ptr = outPtrTmp + outIncZ;
              *pixel += *ptr;
              for (int i = 1; i < numComponents; i++)
                {
                pixel[i] += ptr[i];
                }
              }
            if (count > 0)
              {
              hit = true;
              *maskPtrTmp = newMask;
              Convert((*pixel)/count, outPtrTmp);
              for (int i = 1; i < numComponents; i++)
                {
                Convert(pixel[i]/count, &outPtrTmp[i]);
                }
              }
            }
          maskPtrTmp++;
          outPtrTmp += numComponents;
          }
        }
      }
    lastMask = newMask;
    newMask = 3 - newMask;
    // break if no voxels were set
    if (!hit)
      {
      break;
      }
    }

  delete [] pixel;
}

} // end anonymous namespace

//----------------------------------------------------------------------------
int vtkImageSpread::RequestUpdateExtent(
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
int vtkImageSpread::RequestData(
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

  // we need all the voxels that might be connected to the seed
  int extent[6];
  inData->GetExtent(extent);

  // voxels outside the stencil extent can be excluded
  if (stencil)
    {
    int stencilExtent[6];
    stencil->GetExtent(stencilExtent);
    if (!vtkISF::IntersectExtents(extent, stencilExtent, extent))
      {
      // if stencil doesn't overlap the input, return
      return 1;
      }
    }

  // create and clear the image mask
  size_t size = (extent[1] - extent[0] + 1);
  size *= (extent[3] - extent[2] + 1);
  size *= (extent[5] - extent[4] + 1);
  unsigned char *mask = new unsigned char [size];

  // get scalar pointers
  void *inPtr = inData->GetScalarPointerForExtent(extent);
  void *outPtr = outData->GetScalarPointerForExtent(extent);

  int rval = 1;

  switch (inData->GetScalarType())
    {
    vtkTemplateAliasMacro(
      vtkISF::Execute(this, inData, static_cast<VTK_TT *>(inPtr),
        outData, static_cast<VTK_TT *>(outPtr), mask,
        stencil, extent, this->NumberOfIterations));

    default:
      vtkErrorMacro(<< "Execute: Unknown input ScalarType");
      rval = 0;
    }

  delete [] mask;

  return rval;
}

//----------------------------------------------------------------------------
void vtkImageSpread::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "NumberOfIterations: " << this->NumberOfIterations << "\n";

  os << indent << "StencilConnection: "
     << this->GetStencilConnection() << "\n";
}
