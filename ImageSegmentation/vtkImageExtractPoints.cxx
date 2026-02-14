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
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkImageStencilData.h>
#include <vtkImagePointsIterator.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkSmartPointer.h>
#include <vtkArrayDispatch.h>
#include <vtkDataArrayAccessor.h>

vtkStandardNewMacro(vtkImageExtractPoints);

//----------------------------------------------------------------------------
vtkImageExtractPoints::vtkImageExtractPoints()
{
  this->OutputPointsPrecision = 2;

  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkImageExtractPoints::~vtkImageExtractPoints()
{
}

//----------------------------------------------------------------------------
void vtkImageExtractPoints::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "OutputPointsPrecision: "
     << this->OutputPointsPrecision << "\n";
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
  this->SetInputData(1, stencil);
}

//----------------------------------------------------------------------------
int vtkImageExtractPoints::FillInputPortInformation(
  int port, vtkInformation *info)
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
int vtkImageExtractPoints::FillOutputPortInformation(
  int port, vtkInformation* info)
{
  if (port == 0)
  {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  }

  return 1;
}

//----------------------------------------------------------------------------
int vtkImageExtractPoints::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *vtkNotUsed(outputVector))
{
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
vtkIdType vtkImageExtractPointsCount(
  vtkImageData *inData, vtkImageStencilData *stencil, const int extent[6])
{
  vtkIdType count = 0;
  vtkImageRegionIteratorBase inIter(inData, extent, stencil);
  for (; !inIter.IsAtEnd(); inIter.NextSpan())
  {
    if (inIter.IsInStencil())
    {
      count += inIter.SpanEndId() - inIter.GetId();
    }
  }
  return count;
}

//----------------------------------------------------------------------------
// Worker for vtkArrayDispatch
struct ExtractPointsWorker
{
  vtkImageExtractPoints* Self;
  vtkImageData* InData;
  const int* Extent;
  vtkImageStencilData* Stencil;
  vtkPointData* OutPD;

  template <typename ArrayType>
  void operator()(ArrayType* pointsArray)
  {
    // Use vtk::GetAPIType to get float or double from the array
    using T = vtk::GetAPIType<ArrayType>;
    T* outPoints = pointsArray->GetPointer(0);

    vtkDataArray* inScalars = this->InData->GetPointData()->GetScalars();
    vtkDataArray* outScalars = this->OutPD->GetScalars();
    vtkImagePointsIterator inIter(this->InData, this->Extent, this->Stencil, this->Self, 0);
    
    vtkIdType outId = 0;
    int pixelBytes = inScalars->GetNumberOfComponents() * inScalars->GetDataTypeSize();

    while (!inIter.IsAtEnd())
    {
      if (inIter.IsInStencil())
      {
        vtkIdType n = inIter.SpanEndId() - inIter.GetId();
        void* inPtr = inIter.GetVoidPointer(inScalars, inIter.GetId());
        void* outPtr = outScalars->GetVoidPointer(outId);
        
        memcpy(outPtr, inPtr, n * pixelBytes);
        outId += n;

        for (vtkIdType i = 0; i < n; i++)
        {
          inIter.GetPosition(outPoints);
          outPoints += 3;
          inIter.Next();
        }
      }
      else
      {
        inIter.NextSpan();
      }
    }
  }
};

} // end anonymous namespace

//----------------------------------------------------------------------------
int vtkImageExtractPoints::RequestData(
  vtkInformation*,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  vtkInformation* info = inputVector[0]->GetInformationObject(0);
  vtkInformation *stencilInfo = inputVector[1]->GetInformationObject(0);
  vtkImageData *inData = vtkImageData::SafeDownCast(
    info->Get(vtkDataObject::DATA_OBJECT()));

  vtkImageStencilData* stencil = nullptr;
  if (stencilInfo)
  {
    stencil = static_cast<vtkImageStencilData *>(
      stencilInfo->Get(vtkDataObject::DATA_OBJECT()));
  }

  int pointsType = VTK_DOUBLE;
  if (this->OutputPointsPrecision == 0)
  {
    pointsType = VTK_FLOAT;
  }

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkPolyData *outData = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  const int *extent = inData->GetExtent();
  vtkIdType numPoints = vtkImageExtractPointsCount(inData, stencil, extent);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetDataType(pointsType);
  points->SetNumberOfPoints(numPoints);
  outData->SetPoints(points);

  vtkDataArray *inScalars = inData->GetPointData()->GetScalars();
  vtkPointData *outPD = outData->GetPointData();
  
  vtkSmartPointer<vtkDataArray> scalars;
  scalars.TakeReference(vtkDataArray::CreateDataArray(inScalars->GetDataType()));
  scalars->SetNumberOfComponents(inScalars->GetNumberOfComponents());
  scalars->SetNumberOfTuples(numPoints);
  outPD->SetScalars(scalars);

  // Set up the worker
  ExtractPointsWorker worker;
  worker.Self = this;
  worker.InData = inData;
  worker.Extent = extent;
  worker.Stencil = stencil;
  worker.OutPD = outPD;

  // Use DispatchByValueType to handle float/double safely
  using Dispatcher = vtkArrayDispatch::DispatchByValueType<vtkArrayDispatch::Reals>;
  vtkDataArray* pointsData = points->GetData();

  if (!Dispatcher::Execute(pointsData, worker))
  {
    // Explicitly handle fallback to concrete types to avoid vtkDataArray abstract error
    if (auto* fa = vtkFloatArray::SafeDownCast(pointsData))
    {
      worker(fa);
    }
    else if (auto* da = vtkDoubleArray::SafeDownCast(pointsData))
    {
      worker(da);
    }
    else
    {
      vtkErrorMacro("Unsupported points precision type.");
      return 0;
    }
  }

  return 1;
}
