/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageCorrelationRatio.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageCorrelationRatio.h"

#include "vtkObjectFactory.h"
#include "vtkMath.h"
#include "vtkImageData.h"
#include "vtkImageStencilData.h"
#include "vtkImageStencilIterator.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkMultiThreader.h"
#include "vtkTemplateAliasMacro.h"
#include "vtkPointData.h"
#include "vtkVersion.h"

#include <math.h>

vtkStandardNewMacro(vtkImageCorrelationRatio);

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageCorrelationRatio::vtkImageCorrelationRatio()
{
  this->DataRange[0] = 0;
  this->DataRange[1] = 255;

  this->NumberOfBins = 256;
  this->BinOrigin = 0.0;
  this->BinSpacing = 1.0;

  this->CorrelationRatio = 0.0;

  this->SetNumberOfInputPorts(3);
  this->SetNumberOfOutputPorts(0);
}

//----------------------------------------------------------------------------
vtkImageCorrelationRatio::~vtkImageCorrelationRatio()
{
}

//----------------------------------------------------------------------------
void vtkImageCorrelationRatio::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Stencil: " << this->GetStencil() << "\n";

  os << indent << "DataRange: " << this->DataRange[0] << " "
     << this->DataRange[1] << "\n";

  os << indent << "CorrelationRatio: " << this->CorrelationRatio << "\n";
}

//----------------------------------------------------------------------------
void vtkImageCorrelationRatio::SetStencilData(vtkImageStencilData *stencil)
{
#if VTK_MAJOR_VERSION >= 6
  this->SetInputData(2, stencil);
#else
  this->SetInput(2, stencil);
#endif
}

//----------------------------------------------------------------------------
vtkImageStencilData *vtkImageCorrelationRatio::GetStencil()
{
  if (this->GetNumberOfInputConnections(2) < 1)
    {
    return NULL;
    }
  return vtkImageStencilData::SafeDownCast(
    this->GetExecutive()->GetInputData(2, 0));
}

//----------------------------------------------------------------------------
int vtkImageCorrelationRatio::FillInputPortInformation(
  int port, vtkInformation *info)
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
int vtkImageCorrelationRatio::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* vtkNotUsed(info))
{
  return 1;
}

//----------------------------------------------------------------------------
int vtkImageCorrelationRatio::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *vtkNotUsed(outputVector))
{
  return 1;
}

//----------------------------------------------------------------------------
int vtkImageCorrelationRatio::RequestUpdateExtent(
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
template<class T1, class T2, class T3>
void vtkImageCorrelationRatioExecute(
  vtkImageCorrelationRatio *self,
  vtkImageData *inData0, vtkImageData *inData1, vtkImageStencilData *stencil,
  T1 *inPtr, T2 *inPtr1, int extent[6],
  T3 *outPtr, int numBins, double binOrigin, double binSpacing,
  int threadId)
{
  vtkImageStencilIterator<T1>
    inIter(inData0, stencil, extent, ((threadId == 0) ? self : NULL));
  vtkImageStencilIterator<T2>
    inIter1(inData1, stencil, extent, NULL);

  int pixelInc = inData0->GetNumberOfScalarComponents();
  int pixelInc1 = inData1->GetNumberOfScalarComponents();

  double xmin = 0;
  double xmax = numBins - 1;
  double xshift = -binOrigin;
  double xscale = 1.0/binSpacing;

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
        double x = *inPtr;

        x += xshift;
        x *= xscale;

        x = (x > xmin ? x : xmin);
        x = (x < xmax ? x : xmax);

        int xi = static_cast<int>(x + 0.5);
        T3 *outPtr1 = outPtr + 3*xi;
        T3 y = *inPtr1;
        outPtr1[0]++;
        outPtr1[1] += y;
        outPtr1[2] += y*y;

        inPtr += pixelInc;
        inPtr1 += pixelInc1;
        }
      }
    inIter.NextSpan();
    inIter1.NextSpan();
    }
}

//----------------------------------------------------------------------------
template<class T1, class T2, class T3>
void vtkImageCorrelationRatioExecuteInt(
  vtkImageCorrelationRatio *self,
  vtkImageData *inData0, vtkImageData *inData1, vtkImageStencilData *stencil,
  T1 *inPtr, T2 *inPtr1, int extent[6],
  T3 *outPtr, int numBins, int binOrigin, int binSpacing,
  int threadId)
{
  vtkImageStencilIterator<T1>
    inIter(inData0, stencil, extent, ((threadId == 0) ? self : NULL));
  vtkImageStencilIterator<T2>
    inIter1(inData1, stencil, extent, NULL);

  int pixelInc = inData0->GetNumberOfScalarComponents();
  int pixelInc1 = inData1->GetNumberOfScalarComponents();

  int xmin = 0;
  int xmax = numBins - 1;

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
        int xi = *inPtr;

        xi -= binOrigin;
        xi /= binSpacing;
        xi = (xi > xmin ? xi : xmin);
        xi = (xi < xmax ? xi : xmax);

        T3 *outPtr1 = outPtr + 3*xi;
        T3 y = *inPtr1;
        outPtr1[0]++;
        outPtr1[1] += y;
        outPtr1[2] += y*y;

        inPtr += pixelInc;
        inPtr1 += pixelInc1;
        }
      }
    inIter.NextSpan();
    inIter1.NextSpan();
    }
}

} // end anonymous namespace

//----------------------------------------------------------------------------
// override from vtkThreadedImageAlgorithm to customize the multithreading
int vtkImageCorrelationRatio::RequestData(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // set the scalar information
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *inScalarInfo = vtkDataObject::GetActiveFieldInformation(
    inInfo, vtkDataObject::FIELD_ASSOCIATION_POINTS,
    vtkDataSetAttributes::SCALARS);
  int scalarType = inScalarInfo->Get(vtkDataObject::FIELD_ARRAY_TYPE());

  // compute the array size for the partial sums
  if (scalarType == VTK_DOUBLE || scalarType == VTK_FLOAT)
    {
    this->NumberOfBins = 4096;
    this->BinOrigin = this->DataRange[0];
    double l = this->DataRange[1] - this->DataRange[0];
    this->BinSpacing = l / (this->NumberOfBins - 1);
    }
  else
    {
    this->NumberOfBins = 4096;
    this->BinOrigin = static_cast<int>(this->DataRange[0]);
    int l = static_cast<int>(this->DataRange[1]) - this->BinOrigin;
    if (l < this->NumberOfBins)
      {
      this->NumberOfBins = l + 1;
      }
    this->BinSpacing = (l + this->NumberOfBins)/this->NumberOfBins;
    }

  // specifics for vtkImageCorrelationRatio:
  // allocate workspace for each thread
  vtkIdType memSize = 3*this->NumberOfBins;

  int nThreads = this->GetNumberOfThreads();
  for (int k = 0; k < nThreads; k++)
    {
    this->ThreadOutput[k] = new double[memSize];
    this->ThreadExecuted[k] = false;
    }

  // defer to vtkThreadedImageAlgorithm
  this->Superclass::RequestData(request, inputVector, outputVector);

  // get the dimensions of the joint histogram
  int nx = this->NumberOfBins;

  // piece together the partial sums from each thread
  double n = 0;
  double ySum = 0;
  double yySum = 0;
  double viSum = 0;
  for (int ix = 0; ix < nx; ++ix)
    {
    double ni = 0;
    double yi = 0;
    double yyi = 0;

    for (int j = 0; j < nThreads; j++)
      {
      if (this->ThreadExecuted[j])
        {
        double *outPtr1 = this->ThreadOutput[j] + 3*ix;
        ni += outPtr1[0];
        yi += outPtr1[1];
        yyi += outPtr1[2];
        }
      }

    if (ni > 0)
      {
      // compute the sums for the total variance
      n += ni;
      ySum += yi;
      yySum += yyi;

      // compute the sum for the numerator
      viSum += (yyi - yi*yi/ni);
      }
    }

  // compute the total variance (the denominator)
  double v = 0;
  if (n > 0)
    {
    v = (yySum - ySum*ySum/n);
    }

  // compute the correlation ratio
  double correlationRatio = 0;
  if (v > 0)
    {
    correlationRatio = 1.0 - viSum/v;
    }

  // delete the temporary memory
  for (int j = 0; j < nThreads; j++)
    {
    delete [] this->ThreadOutput[j];
    }

  // output values
  this->CorrelationRatio = correlationRatio;

  return 1;
}

//----------------------------------------------------------------------------
// turn off 64-bit ints when templating over all types
# undef VTK_USE_INT64
# define VTK_USE_INT64 0
# undef VTK_USE_UINT64
# define VTK_USE_UINT64 0

namespace {

//----------------------------------------------------------------------------
// turn off floats in vtkTemplateAliasMacro specifically for this execute
# undef VTK_USE_FLOAT32
# define VTK_USE_FLOAT32 0
# undef VTK_USE_FLOAT64
# define VTK_USE_FLOAT64 0

template<class T1>
void vtkImageCorrelationRatioExecute1Int(
  vtkImageCorrelationRatio *self,
  vtkImageData *inData0, vtkImageData *inData1, vtkImageStencilData *stencil,
  void *inPtr, T1 *inPtr1, int extent[6],
  double *outPtr, int numBins, int binOrigin, int binSpacing,
  int threadId)
{
  switch (inData0->GetScalarType())
    {
    vtkTemplateAliasMacro(
      vtkImageCorrelationRatioExecuteInt(
        self, inData0, inData1, stencil,
        static_cast<VTK_TT *>(inPtr), inPtr1, extent,
        outPtr, numBins, binOrigin, binSpacing, threadId));
    default:
      vtkErrorWithObjectMacro(self, "Execute: Unknown input ScalarType");
    }
}

// turn floats back on
# undef VTK_USE_FLOAT32
# define VTK_USE_FLOAT32 1
# undef VTK_USE_FLOAT64
# define VTK_USE_FLOAT64 1

//----------------------------------------------------------------------------
template<class T1>
void vtkImageCorrelationRatioExecute1(
  vtkImageCorrelationRatio *self,
  vtkImageData *inData0, vtkImageData *inData1, vtkImageStencilData *stencil,
  void *inPtr, T1 *inPtr1, int extent[6],
  double *outPtr, int numBins, double binOrigin, double binSpacing,
  int threadId)
{
  if (inData0->GetScalarType() == VTK_FLOAT)
    {
    vtkImageCorrelationRatioExecute(
      self, inData0, inData1, stencil,
      static_cast<float *>(inPtr), inPtr1, extent,
      outPtr, numBins, binOrigin, binSpacing, threadId);
    }
  else if (inData0->GetScalarType() == VTK_DOUBLE)
    {
    vtkImageCorrelationRatioExecute(
      self, inData0, inData1, stencil,
      static_cast<double *>(inPtr), inPtr1, extent,
      outPtr, numBins, binOrigin, binSpacing, threadId);
    }
}

} // end anonymous namespace

//----------------------------------------------------------------------------
// This method is passed a input and output region, and executes the filter
// algorithm to fill the output from the input.
// It just executes a switch statement to call the correct function for
// the regions data types.
void vtkImageCorrelationRatio::ThreadedRequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector),
  vtkImageData ***vtkNotUsed(inData),
  vtkImageData **vtkNotUsed(outData),
  int extent[6], int threadId)
{
  this->ThreadExecuted[threadId] = true;
  double *outPtr = this->ThreadOutput[threadId];

  // initialize the partial sums to zero
  vtkIdType outCount = 3*this->NumberOfBins;
  double *outPtr1 = outPtr;
  do { *outPtr1++ = 0; } while (--outCount > 0);

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

  double binOrigin = this->BinOrigin;
  double binSpacing = this->BinSpacing;
  int numBins = this->NumberOfBins;

  if (inData0->GetScalarType() != VTK_FLOAT &&
      inData0->GetScalarType() != VTK_DOUBLE)
    {
    switch (inData1->GetScalarType())
      {
      vtkTemplateAliasMacro(
        vtkImageCorrelationRatioExecute1Int(
          this, inData0, inData1, stencil,
          inPtr0, static_cast<VTK_TT *>(inPtr1),
          extent, outPtr, numBins,
          static_cast<int>(binOrigin), static_cast<int>(binSpacing),
          threadId));
      default:
        vtkErrorMacro(<< "Execute: Unknown ScalarType");
      }
    }
  else
    {
    switch (inData1->GetScalarType())
      {
      vtkTemplateAliasMacro(
        vtkImageCorrelationRatioExecute1(
          this, inData0, inData1, stencil,
          inPtr0, static_cast<VTK_TT *>(inPtr1),
          extent, outPtr, numBins, binOrigin, binSpacing,
          threadId));
      default:
        vtkErrorMacro(<< "Execute: Unknown ScalarType");
      }
    }
}
