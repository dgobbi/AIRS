/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageCrossCorrelation.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageCrossCorrelation.h"

#include "vtkObjectFactory.h"
#include "vtkImageData.h"
#include "vtkImageStencilData.h"
#include "vtkImageStencilIterator.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkMultiThreader.h"
#include "vtkTemplateAliasMacro.h"

// turn off 64-bit ints when templating over all types
# undef VTK_USE_INT64
# define VTK_USE_INT64 0
# undef VTK_USE_UINT64
# define VTK_USE_UINT64 0

#include <math.h>

vtkStandardNewMacro(vtkImageCrossCorrelation);

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageCrossCorrelation::vtkImageCrossCorrelation()
{
  this->CrossCorrelation = 0.0;
  this->NormalizedCrossCorrelation = 0.0;

  this->SetNumberOfInputPorts(3);
  this->SetNumberOfOutputPorts(0);
}

//----------------------------------------------------------------------------
vtkImageCrossCorrelation::~vtkImageCrossCorrelation()
{
}

//----------------------------------------------------------------------------
void vtkImageCrossCorrelation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Stencil: " << this->GetStencil() << "\n";

  os << indent << "CrossCorrelation: " << this->CrossCorrelation << "\n";
  os << indent << "NormalizedCrossCorrelation: "
     << this->NormalizedCrossCorrelation << "\n";
}

//----------------------------------------------------------------------------
void vtkImageCrossCorrelation::SetStencil(vtkImageStencilData *stencil)
{
  this->SetInput(2, stencil);
}

//----------------------------------------------------------------------------
vtkImageStencilData *vtkImageCrossCorrelation::GetStencil()
{
  if (this->GetNumberOfInputConnections(2) < 1)
    {
    return NULL;
    }
  return vtkImageStencilData::SafeDownCast(
    this->GetExecutive()->GetInputData(2, 0));
}

//----------------------------------------------------------------------------
int vtkImageCrossCorrelation::FillInputPortInformation(int port, vtkInformation *info)
{
  if (port == 0 || port == 1)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    }
  else if (port == 2)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageStencilData");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }

  return 1;
}

//----------------------------------------------------------------------------
int vtkImageCrossCorrelation::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* vtkNotUsed(info))
{
  return 1;
}

//----------------------------------------------------------------------------
int vtkImageCrossCorrelation::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *vtkNotUsed(outputVector))
{
  return 1;
}

//----------------------------------------------------------------------------
int vtkImageCrossCorrelation::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector))
{
  int inExt0[6], inExt1[6];
  vtkInformation *inInfo0 = inputVector[0]->GetInformationObject(0);
  vtkInformation *inInfo1 = inputVector[1]->GetInformationObject(0);

  inInfo0->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), inExt0);
  inInfo1->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), inExt1);

  inInfo0->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inExt0, 6);
  inInfo1->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inExt1, 6);

  // need to set the stencil update extent to the input extent
  if (this->GetNumberOfInputConnections(2) > 0)
    {
    vtkInformation *stencilInfo = inputVector[2]->GetInformationObject(0);
    stencilInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
                     inExt0, 6);
    }

  return 1;
}

// begin anonymous namespace
namespace {

//----------------------------------------------------------------------------
template<class T1, class T2>
void vtkImageCrossCorrelationExecute(
  vtkImageCrossCorrelation *self,
  vtkImageData *inData0, vtkImageData *inData1, vtkImageStencilData *stencil,
  T1 *inPtr, T2 *inPtr1, int extent[6], double output[4],
  int threadId)
{
  vtkImageStencilIterator<T1>
    inIter(inData0, stencil, extent, ((threadId == 0) ? self : NULL));
  vtkImageStencilIterator<T2>
    inIter1(inData1, stencil, extent, NULL);

  double xSum = 0;
  double ySum = 0;
  double xySum = 0;
  double count = 0;

  // iterate over all spans in the stencil
  while (!inIter.IsAtEnd())
    {
    if (inIter.IsInStencil())
      {
      inPtr = inIter.BeginSpan();
      T1 *inPtrEnd = inIter.EndSpan();
      inPtr1 = inIter1.BeginSpan();

      // iterate over all voxels in the span
      while (inPtr != inPtrEnd)
        {
        double x = *inPtr++;
        double y = *inPtr1++;

        xSum += x*x;
        ySum += y*y;
        xySum += x*y;
        count++;
        }
      }
    inIter.NextSpan();
    inIter1.NextSpan();
    }

  output[0] = xSum;
  output[1] = ySum;
  output[2] = xySum;
  output[3] = count;
}

//----------------------------------------------------------------------------
template<class T1>
void vtkImageCrossCorrelationExecute1(
  vtkImageCrossCorrelation *self,
  vtkImageData *inData0, vtkImageData *inData1, vtkImageStencilData *stencil,
  T1 *inPtr, void *inPtr1, int extent[6], double output[4], int threadId)
{
  switch (inData1->GetScalarType())
    {
    vtkTemplateAliasMacro(
      vtkImageCrossCorrelationExecute(
        self, inData0, inData1, stencil,
        inPtr, static_cast<VTK_TT *>(inPtr1), extent, output, threadId));
    default:
      vtkErrorWithObjectMacro(self, "Execute: Unknown input ScalarType");
    }
}

} // end anonymous namespace

//----------------------------------------------------------------------------
// override from vtkThreadedImageAlgorithm to customize the multithreading
int vtkImageCrossCorrelation::RequestData(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // specifics for vtkImageCrossCorrelation:
  // allocate workspace for each thread

  int n = this->GetNumberOfThreads();
  for (int k = 0; k < n; k++)
    {
    this->ThreadOutput[k][0] = 0;
    this->ThreadOutput[k][1] = 0;
    this->ThreadOutput[k][2] = 0;
    this->ThreadOutput[k][3] = 0;
    }

  // defer to vtkThreadedImageAlgorithm
  this->Superclass::RequestData(request, inputVector, outputVector);

  // various variables for computing the mutual information
  double xSum = 0;
  double ySum = 0;
  double xySum = 0;
  double count = 0;

  // add the contribution from thread j
  for (int j = 0; j < n; j++)
    {
    xSum += this->ThreadOutput[j][0];
    ySum += this->ThreadOutput[j][1];
    xySum += this->ThreadOutput[j][2];
    count += this->ThreadOutput[j][3];
    }

  // minimum possible values
  double crossCorrelation = 0;
  double normalizedCrossCorrelation = 1.0;

  if (count > 0)
    {
    crossCorrelation = xySum;

    if (xSum > 0 && ySum > 0)
      {
      normalizedCrossCorrelation = xySum/sqrt(xSum*ySum);
      }
    }

  // output values
  this->CrossCorrelation = crossCorrelation;
  this->NormalizedCrossCorrelation = normalizedCrossCorrelation;

  return 1;
}

//----------------------------------------------------------------------------
// This method is passed a input and output region, and executes the filter
// algorithm to fill the output from the input.
// It just executes a switch statement to call the correct function for
// the regions data types.
void vtkImageCrossCorrelation::ThreadedRequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector),
  vtkImageData ***vtkNotUsed(inData),
  vtkImageData **vtkNotUsed(outData),
  int extent[6], int threadId)
{
  vtkInformation *inInfo0 = inputVector[0]->GetInformationObject(0);
  vtkInformation *inInfo1 = inputVector[1]->GetInformationObject(0);

  vtkImageData *inData0 = vtkImageData::SafeDownCast(
    inInfo0->Get(vtkDataObject::DATA_OBJECT()));
  vtkImageData *inData1 = vtkImageData::SafeDownCast(
    inInfo1->Get(vtkDataObject::DATA_OBJECT()));

  // make sure execute extent is not beyond the extent of any input
  int inExt0[6], inExt1[6];
  inData0->GetExtent(inExt0);
  inData1->GetExtent(inExt1);

  for (int i = 0; i < 6; i += 2)
    {
    int j = i + 1;
    extent[i] = ((extent[i] > inExt0[i]) ? extent[i] : inExt0[i]);
    extent[i] = ((extent[i] > inExt1[i]) ? extent[i] : inExt1[i]);
    extent[j] = ((extent[j] < inExt0[j]) ? extent[j] : inExt0[j]);
    extent[j] = ((extent[j] < inExt1[j]) ? extent[j] : inExt1[j]);
    if (extent[i] > extent[j])
      {
      return;
      }
    }

  void *inPtr0 = inData0->GetScalarPointerForExtent(extent);
  void *inPtr1 = inData1->GetScalarPointerForExtent(extent);

  vtkImageStencilData *stencil = this->GetStencil();

  switch (inData0->GetScalarType())
    {
    vtkTemplateAliasMacro(
      vtkImageCrossCorrelationExecute1(
        this, inData0, inData1, stencil,
        static_cast<VTK_TT *>(inPtr0), inPtr1, extent,
        this->ThreadOutput[threadId], threadId));
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
    }
}
