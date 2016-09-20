/*=========================================================================

  Module: vtkImageCrossCorrelation.cxx

  Copyright (c) 2006 Atamai, Inc.
  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageCrossCorrelation.h"

#include "vtkImageSimilarityMetricInternals.h"

#include <vtkObjectFactory.h>
#include <vtkImageData.h>
#include <vtkImageStencilData.h>
#include <vtkImageStencilIterator.h>
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

vtkStandardNewMacro(vtkImageCrossCorrelation);

//----------------------------------------------------------------------------
// Data needed for each thread.
class vtkImageCrossCorrelationThreadData
{
public:
  vtkImageCrossCorrelationThreadData()
  {
    Data[0] = Data[1] = Data[2] = Data[3] = Data[4] = Data[5] = 0.0;
  }

  double Data[6];
};

class vtkImageCrossCorrelationTLS
  : public vtkImageSimilarityMetricTLS<vtkImageCrossCorrelationThreadData>
{
};

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageCrossCorrelation::vtkImageCrossCorrelation()
{
  this->Metric = NCC;
  this->CrossCorrelation = 0.0;
  this->NormalizedCrossCorrelation = 0.0;

  this->ThreadData = 0;
}

//----------------------------------------------------------------------------
vtkImageCrossCorrelation::~vtkImageCrossCorrelation()
{
}

//----------------------------------------------------------------------------
void vtkImageCrossCorrelation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "CrossCorrelation: " << this->CrossCorrelation << "\n";
  os << indent << "NormalizedCrossCorrelation: "
     << this->NormalizedCrossCorrelation << "\n";
  os << indent << "Metric: "
     << (this->Metric == NCC ? "NormalizedCrossCorrelation\n" :
                               "CrossCorrelation\n");
}

// begin anonymous namespace
namespace {

//----------------------------------------------------------------------------
template<class T1, class T2>
void vtkImageCrossCorrelationExecute(
  vtkImageCrossCorrelation *self,
  vtkImageData *inData0, vtkImageData *inData1, vtkImageStencilData *stencil,
  T1 *inPtr, T2 *inPtr1, const int extent[6], double output[6],
  vtkIdType pieceId)
{
  int *ext = const_cast<int *>(extent);
  vtkImageStencilIterator<T1>
    inIter(inData0, stencil, ext, ((pieceId == 0) ? self : NULL));
  vtkImageStencilIterator<T2>
    inIter1(inData1, stencil, ext, NULL);

  int pixelInc = inData0->GetNumberOfScalarComponents();
  int pixelInc1 = inData1->GetNumberOfScalarComponents();

  double xSum = 0;
  double ySum = 0;
  double xxSum = 0;
  double yySum = 0;
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
        double x = *inPtr;
        double y = *inPtr1;

        xSum += x;
        ySum += y;
        xxSum += x*x;
        yySum += y*y;
        xySum += x*y;
        count++;

        inPtr += pixelInc;
        inPtr1 += pixelInc1;
        }
      }
    inIter.NextSpan();
    inIter1.NextSpan();
    }

  output[0] = xSum;
  output[1] = ySum;
  output[2] = xxSum;
  output[3] = yySum;
  output[4] = xySum;
  output[5] = count;
}

//----------------------------------------------------------------------------
template<class T1>
void vtkImageCrossCorrelationExecute1(
  vtkImageCrossCorrelation *self,
  vtkImageData *inData0, vtkImageData *inData1,
  vtkImageStencilData *stencil, T1 *inPtr, void *inPtr1,
  const int extent[6], double output[6], vtkIdType pieceId)
{
  switch (inData1->GetScalarType())
    {
    vtkTemplateAliasMacro(
      vtkImageCrossCorrelationExecute(
        self, inData0, inData1, stencil,
        inPtr, static_cast<VTK_TT *>(inPtr1), extent, output, pieceId));
    default:
      vtkErrorWithObjectMacro(self, "Execute: Unknown input ScalarType");
    }
}

} // end anonymous namespace

//----------------------------------------------------------------------------
int vtkImageCrossCorrelation::RequestData(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // create the thread-local object
  vtkImageCrossCorrelationTLS tlocal;
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
void vtkImageCrossCorrelation::PieceRequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector),
  const int extent[6], vtkIdType pieceId)
{
  vtkImageCrossCorrelationThreadData *threadLocal =
    &this->ThreadData->Local(pieceId);

  double *outPtr = threadLocal->Data;

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

  switch (inData0->GetScalarType())
    {
    vtkTemplateAliasMacro(
      vtkImageCrossCorrelationExecute1(
        this, inData0, inData1, stencil,
        static_cast<VTK_TT *>(inPtr0), inPtr1, extent,
        outPtr, pieceId));
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
    }
}

//----------------------------------------------------------------------------
void vtkImageCrossCorrelation::ReduceRequestData(
  vtkInformation *, vtkInformationVector **, vtkInformationVector *)
{
  // various variables for computing the mutual information
  double xSum = 0.0;
  double ySum = 0.0;
  double xxSum = 0.0;
  double yySum = 0.0;
  double xySum = 0.0;
  double count = 0.0;

  // add the contributions from all threads
  for (vtkImageCrossCorrelationTLS::iterator
       iter = this->ThreadData->begin();
       iter != this->ThreadData->end(); ++iter)
    {
    double *data = iter->Data;
    xSum += data[0];
    ySum += data[1];
    xxSum += data[2];
    yySum += data[3];
    xySum += data[4];
    count += data[5];
    }

  // minimum possible values
  double crossCorrelation = 0.0;
  double normalizedCrossCorrelation = 1.0;

  if (count > 0)
    {
    crossCorrelation = (xySum - xSum*ySum/count)/count;

    if (xxSum > 0 && yySum > 0)
      {
      normalizedCrossCorrelation = (xySum - xSum*ySum/count)/
        sqrt((xxSum - xSum*xSum/count)*(yySum - ySum*ySum/count));

      // was double precision able to capture the sums exactly?
      if (xxSum > 1e16 || yySum > 1e16)
        {
        vtkWarningMacro("Possible incorrect result due to "
                        "insufficient precision in subtraction");
        }
      }
    }

  // output values
  this->CrossCorrelation = crossCorrelation;
  this->NormalizedCrossCorrelation = normalizedCrossCorrelation;

  if (this->Metric == NCC)
    {
    this->SetCost(-normalizedCrossCorrelation);
    }
  else
    {
    this->SetCost(-crossCorrelation);
    }
}
