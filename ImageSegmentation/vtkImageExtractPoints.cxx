/*=========================================================================

  Copyright (c) 2015,2016 David Gobbi
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

  * Neither the name of David Gobbi nor the names of any contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#include "vtkImageExtractPoints.h"

#include <vtkObjectFactory.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkIdTypeArray.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkImageStencilData.h>
#include <vtkImageRegionIterator.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTemplateAliasMacro.h>

// turn off 64-bit ints when templating over all types
# undef VTK_USE_INT64
# define VTK_USE_INT64 0
# undef VTK_USE_UINT64
# define VTK_USE_UINT64 0

vtkStandardNewMacro(vtkImageExtractPoints);

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageExtractPoints::vtkImageExtractPoints()
{
  this->ActiveComponent = -1;
  this->Points = vtkPoints::New();
  this->Points->SetDataTypeToDouble();
  this->PointIds = vtkIdTypeArray::New();
  this->Scalars = vtkDoubleArray::New();

  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkImageExtractPoints::~vtkImageExtractPoints()
{
  if (this->Points)
    {
    this->Points->Delete();
    }
  if (this->PointIds)
    {
    this->PointIds->Delete();
    }
  if (this->Scalars)
    {
    this->Scalars->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkImageExtractPoints::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "ActiveComponent: " << this->ActiveComponent << "\n";
  os << indent << "Points: " << this->Points << "\n";
  os << indent << "PointIds: " << this->PointIds << "\n";
  os << indent << "Scalars: " << this->Scalars << "\n";
}

//----------------------------------------------------------------------------
vtkPoints *vtkImageExtractPoints::GetPoints()
{
  return this->Points;
}

//----------------------------------------------------------------------------
vtkIdTypeArray *vtkImageExtractPoints::GetPointIds()
{
  return this->PointIds;
}

//----------------------------------------------------------------------------
vtkDataArray *vtkImageExtractPoints::GetScalars()
{
  return this->Scalars;
}

//----------------------------------------------------------------------------
void vtkImageExtractPoints::SetStencilConnection(
  vtkAlgorithmOutput *stencil)
{
  this->SetInputConnection(1, stencil);
}

//----------------------------------------------------------------------------
vtkAlgorithmOutput *vtkImageExtractPoints::GetStencilConnection()
{
  return this->GetInputConnection(1, 0);
}

//----------------------------------------------------------------------------
void vtkImageExtractPoints::SetStencilData(vtkImageStencilData *stencil)
{
#if VTK_MAJOR_VERSION >= 6
  this->SetInputData(1, stencil);
#else
  this->SetInput(1, stencil);
#endif
}

//----------------------------------------------------------------------------
int vtkImageExtractPoints::FillInputPortInformation(int port, vtkInformation *info)
{
  if (port == 0)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    }
  else if (port == 1)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageStencilData");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }

  return 1;
}

//----------------------------------------------------------------------------
int vtkImageExtractPoints::FillOutputPortInformation(int port, vtkInformation* info)
{
  if (port == 0)
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    }

  return 1;
}

//----------------------------------------------------------------------------
int vtkImageExtractPoints::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  if (this->GetNumberOfOutputPorts() > 0)
    {
    int outWholeExt[6] = { 0, 0, 0, 0, 0, 0 };
    double outOrigin[3] = { 0.0, 0.0, 0.0 };
    double outSpacing[3] = { 1.0, 1.0, 1.0 };

    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
                 outWholeExt, 6);

    outInfo->Set(vtkDataObject::ORIGIN(), outOrigin, 3);
    outInfo->Set(vtkDataObject::SPACING(), outSpacing, 3);

    vtkDataObject::SetPointDataActiveScalarInfo(
      outInfo, VTK_UNSIGNED_CHAR, 1);
    }

  return 1;
}

//----------------------------------------------------------------------------
int vtkImageExtractPoints::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector))
{
  int inExt[6];
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

  inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), inExt);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inExt, 6);

  // need to set the stencil update extent to the input extent
  if (this->GetNumberOfInputConnections(1) > 0)
    {
    vtkInformation *stencilInfo = inputVector[1]->GetInformationObject(0);
    stencilInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
                     inExt, 6);
    }

  return 1;
}

namespace {

//----------------------------------------------------------------------------
vtkIdType vtkImageExtractPointsCount(vtkImageStencilData *stencil, int extent[6])
{
  // count the number of voxels inside the stencil
  vtkIdType count = 0;
  for (int k = extent[4]; k <= extent[5]; k++)
    {
    for (int j = extent[2]; j <= extent[3]; j++)
      {
      int iter = 0;
      int r1, r2;
      while (stencil->GetNextExtent(r1, r2, extent[0], extent[1], j, k, iter))
        {
        count += r2 - r1 + 1;
        }
      }
    }

  return count;
}

//----------------------------------------------------------------------------
template<class T>
void vtkImageExtractPointsExecute(
  vtkImageData *inData, vtkImageStencilData *stencil,
  const T *inPtr, const int extent[6],
  const double origin[3], const double spacing[3],
  double *outPts, vtkIdType *outPtIds, double *outPtr, int component)
{
  vtkImageRegionIterator<T> inIter(inData, stencil, extent);

  // set up components
  int nc = inData->GetNumberOfScalarComponents();
  int c = component;
  int inInc = nc;
  if (c < 0)
    {
    c = 0;
    inInc = 1;
    }

  // iterate over all spans in the stencil
  while (!inIter.IsAtEnd())
    {
    if (inIter.IsInStencil())
      {
      inPtr = inIter.BeginSpan();
      T *inPtrEnd = inIter.EndSpan();
      vtkIdType ptId = inIter.GetPointId();
      int i = inIter.GetIndexX();
      double y = origin[1] + spacing[1]*inIter.GetIndexY();
      double z = origin[2] + spacing[2]*inIter.GetIndexZ();
      while (inPtr != inPtrEnd)
        {
        outPts[0] = origin[0] + spacing[0]*i;
        outPts[1] = y;
        outPts[2] = z;
        i++;
        *outPtIds = ptId++;
        *outPtr = inPtr[c];
        inPtr += inInc;
        outPts += 3;
        outPtIds++;
        outPtr++;
        }
      }
    inIter.NextSpan();
    }
}

}

//----------------------------------------------------------------------------
int vtkImageExtractPoints::RequestData(
  vtkInformation*,
  vtkInformationVector** inputVector,
  vtkInformationVector*)
{
  vtkInformation* info = inputVector[0]->GetInformationObject(0);
  vtkInformation *stencilInfo = inputVector[1]->GetInformationObject(0);
  vtkImageData *inData = vtkImageData::SafeDownCast(
    info->Get(vtkDataObject::DATA_OBJECT()));

  double origin[3], spacing[3];
  inData->GetOrigin(origin);
  inData->GetSpacing(spacing);

  int extent[6];
  inData->GetExtent(extent);
  void *inPtr = inData->GetScalarPointerForExtent(extent);

  vtkImageStencilData* stencil = 0;
  if (stencilInfo)
    {
    stencil = static_cast<vtkImageStencilData *>(
      stencilInfo->Get(vtkDataObject::DATA_OBJECT()));
    }

  this->Points->Initialize();
  this->PointIds->Initialize();
  this->Scalars->Initialize();

  int scalarType = inData->GetScalarType();
  int numComponents = inData->GetNumberOfScalarComponents();
  int component = this->ActiveComponent;

  int outComponents = 1;
  if (component < 0)
    {
    outComponents = numComponents;
    }
  this->Scalars->SetNumberOfComponents(outComponents);

  vtkIdType numVoxels = vtkImageExtractPointsCount(stencil, extent);
  this->Points->SetNumberOfPoints(numVoxels);
  this->PointIds->SetNumberOfTuples(numVoxels);
  this->Scalars->SetNumberOfTuples(numVoxels);

  double *outPtr = this->Scalars->WritePointer(0, numVoxels*outComponents);
  double *outPts = vtkDoubleArray::SafeDownCast(this->Points->GetData())
    ->WritePointer(0, numVoxels*3);
  vtkIdType *outPtIds = this->PointIds->WritePointer(0, numVoxels);

  switch (scalarType)
    {
    vtkTemplateAliasMacro(
      vtkImageExtractPointsExecute(
        inData, stencil, static_cast<const VTK_TT *>(inPtr),
        extent, origin, spacing, outPts, outPtIds, outPtr, component));
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
    }

  return 1;
}
