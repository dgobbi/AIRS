/*=========================================================================

  Module: vtkImageCorrelationRatio.cxx

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageCorrelationRatio.h"

#include "vtkImageSimilarityMetricInternals.h"

#include <vtkObjectFactory.h>
#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkImageStencilData.h>
#include <vtkImageStencilIterator.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTemplateAliasMacro.h>
#include <vtkPointData.h>
#include <vtkVersion.h>

#include <math.h>

vtkStandardNewMacro(vtkImageCorrelationRatio);

//----------------------------------------------------------------------------
// Data needed for each thread.
class vtkImageCorrelationRatioThreadData
{
public:
  vtkImageCorrelationRatioThreadData() : Data(0) {}

  double *Data;
};

class vtkImageCorrelationRatioTLS
  : public vtkImageSimilarityMetricTLS<vtkImageCorrelationRatioThreadData>
{
};

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

  this->ThreadData = 0;
}

//----------------------------------------------------------------------------
vtkImageCorrelationRatio::~vtkImageCorrelationRatio()
{
}

//----------------------------------------------------------------------------
void vtkImageCorrelationRatio::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "DataRange: " << this->DataRange[0] << " "
     << this->DataRange[1] << "\n";

  os << indent << "CorrelationRatio: " << this->CorrelationRatio << "\n";
}

// begin anonymous namespace
namespace {

//----------------------------------------------------------------------------
template<class T1, class T2, class T3>
void vtkImageCorrelationRatioExecute(
  vtkImageCorrelationRatio *self,
  vtkImageData *inData0, vtkImageData *inData1, vtkImageStencilData *stencil,
  T1 *inPtr, T2 *inPtr1, const int extent[6],
  T3 *outPtr, int numBins, double binOrigin, double binSpacing,
  vtkIdType pieceId)
{
  int *ext = const_cast<int *>(extent);
  vtkImageStencilIterator<T1>
    inIter(inData0, stencil, ext, ((pieceId == 0) ? self : NULL));
  vtkImageStencilIterator<T2>
    inIter1(inData1, stencil, ext, NULL);

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
  T1 *inPtr, T2 *inPtr1, const int extent[6],
  T3 *outPtr, int numBins, int binOrigin, int binSpacing,
  vtkIdType pieceId)
{
  int *ext = const_cast<int *>(extent);
  vtkImageStencilIterator<T1>
    inIter(inData0, stencil, ext, ((pieceId == 0) ? self : NULL));
  vtkImageStencilIterator<T2>
    inIter1(inData1, stencil, ext, NULL);

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
int vtkImageCorrelationRatio::RequestData(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // get the scalar information
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

  // create the thread-local object
  vtkImageCorrelationRatioTLS tlocal;
  tlocal.Initialize(this);
  this->ThreadData = &tlocal;

  this->Superclass::RequestData(request, inputVector, outputVector);

  this->ThreadData = 0;

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
  void *inPtr, T1 *inPtr1, const int extent[6],
  double *outPtr, int numBins, int binOrigin, int binSpacing,
  vtkIdType pieceId)
{
  switch (inData0->GetScalarType())
    {
    vtkTemplateAliasMacro(
      vtkImageCorrelationRatioExecuteInt(
        self, inData0, inData1, stencil,
        static_cast<VTK_TT *>(inPtr), inPtr1, extent,
        outPtr, numBins, binOrigin, binSpacing, pieceId));
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
  void *inPtr, T1 *inPtr1, const int extent[6],
  double *outPtr, int numBins, double binOrigin, double binSpacing,
  vtkIdType pieceId)
{
  if (inData0->GetScalarType() == VTK_FLOAT)
    {
    vtkImageCorrelationRatioExecute(
      self, inData0, inData1, stencil,
      static_cast<float *>(inPtr), inPtr1, extent,
      outPtr, numBins, binOrigin, binSpacing, pieceId);
    }
  else if (inData0->GetScalarType() == VTK_DOUBLE)
    {
    vtkImageCorrelationRatioExecute(
      self, inData0, inData1, stencil,
      static_cast<double *>(inPtr), inPtr1, extent,
      outPtr, numBins, binOrigin, binSpacing, pieceId);
    }
}

} // end anonymous namespace

//----------------------------------------------------------------------------
// This method is passed a input and output region, and executes the filter
// algorithm to fill the output from the input.
// It just executes a switch statement to call the correct function for
// the regions data types.
void vtkImageCorrelationRatio::PieceRequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector),
  const int extent[6], vtkIdType pieceId)
{
  vtkImageCorrelationRatioThreadData *threadLocal =
    &this->ThreadData->Local(pieceId);

  double *outPtr = threadLocal->Data;

  if (outPtr == 0)
    {
    // initialize the partial sums to zero
    vtkIdType outCount = 3*this->NumberOfBins;
    threadLocal->Data = new double[outCount];
    outPtr = threadLocal->Data;

    double *outPtr1 = outPtr;
    do { *outPtr1++ = 0; } while (--outCount > 0);
    }

  vtkInformation *inInfo0 = inputVector[0]->GetInformationObject(0);
  vtkInformation *inInfo1 = inputVector[1]->GetInformationObject(0);

  vtkImageData *inData0 = vtkImageData::SafeDownCast(
    inInfo0->Get(vtkDataObject::DATA_OBJECT()));
  vtkImageData *inData1 = vtkImageData::SafeDownCast(
    inInfo1->Get(vtkDataObject::DATA_OBJECT()));

  int *ext = const_cast<int *>(extent);
  void *inPtr0 = inData0->GetScalarPointerForExtent(ext);
  void *inPtr1 = inData1->GetScalarPointerForExtent(ext);

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
          pieceId));
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
          pieceId));
      default:
        vtkErrorMacro(<< "Execute: Unknown ScalarType");
      }
    }
}

//----------------------------------------------------------------------------
void vtkImageCorrelationRatio::ReduceRequestData(
  vtkInformation *, vtkInformationVector **, vtkInformationVector *)
{
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

    for (vtkImageCorrelationRatioTLS::iterator
         iter = this->ThreadData->begin();
         iter != this->ThreadData->end(); ++iter)
      {
      if (iter->Data)
        {
        double *outPtr1 = iter->Data + 3*ix;
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

  for (vtkImageCorrelationRatioTLS::iterator
       iter = this->ThreadData->begin();
       iter != this->ThreadData->end(); ++iter)
    {
    // delete the temporary memory
    delete [] iter->Data;
    }

  // output values
  this->CorrelationRatio = correlationRatio;

  this->SetCost(-correlationRatio);
}
