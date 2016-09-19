/*=========================================================================

  Module: vtkImageSquaredDifference.cxx

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageSquaredDifference.h"

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

vtkStandardNewMacro(vtkImageSquaredDifference);

//----------------------------------------------------------------------------
// Data needed for each thread.
class vtkImageSquaredDifferenceThreadData
{
public:
  vtkImageSquaredDifferenceThreadData() : SumSquares(0.0), Count(0) {}

  double SumSquares;
  vtkIdType Count;
};

class vtkImageSquaredDifferenceTLS
  : public vtkImageSimilarityMetricTLS<vtkImageSquaredDifferenceThreadData>
{
};

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageSquaredDifference::vtkImageSquaredDifference()
{
  this->SquaredDifference = 0.0;
}

//----------------------------------------------------------------------------
vtkImageSquaredDifference::~vtkImageSquaredDifference()
{
}

//----------------------------------------------------------------------------
void vtkImageSquaredDifference::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "SquaredDifference: " << this->SquaredDifference << "\n";
}

// begin anonymous namespace
namespace {

//----------------------------------------------------------------------------
template<class T1, class T2>
void vtkImageSquaredDifferenceExecute(
  vtkImageSquaredDifference *self,
  vtkImageData *inData0, vtkImageData *inData1, vtkImageStencilData *stencil,
  T1 *inPtr, T2 *inPtr1, const int extent[6], vtkIdType pieceId,
  vtkImageSquaredDifferenceThreadData *output)
{
  int *ext = const_cast<int *>(extent);
  vtkImageStencilIterator<T1>
    inIter(inData0, stencil, ext, ((pieceId == 0) ? self : NULL));
  vtkImageStencilIterator<T2>
    inIter1(inData1, stencil, ext, NULL);

  double sqsum = 0;
  vtkIdType count = 0;

  // iterate over all spans in the stencil
  while (!inIter.IsAtEnd())
    {
    if (inIter.IsInStencil())
      {
      inPtr = inIter.BeginSpan();
      T1 *inPtrEnd = inIter.EndSpan();
      inPtr1 = inIter1.BeginSpan();

      double s = 0;

      // iterate over all voxels in the span
      while (inPtr != inPtrEnd)
        {
        double x = *inPtr++;
        double y = *inPtr1++;
        double d = y - x;
        s += d*d;
        }

      count += static_cast<vtkIdType>(inPtrEnd - inPtr);
      sqsum += s;
      }
    inIter.NextSpan();
    inIter1.NextSpan();
    }

  output->SumSquares += sqsum;
  output->Count += count;
}

//----------------------------------------------------------------------------
template<class T1>
void vtkImageSquaredDifferenceExecute1(
  vtkImageSquaredDifference *self,
  vtkImageData *inData0, vtkImageData *inData1, vtkImageStencilData *stencil,
  T1 *inPtr, void *inPtr1, const int extent[6], vtkIdType pieceId,
  vtkImageSquaredDifferenceThreadData *output)
{
  switch (inData1->GetScalarType())
    {
    vtkTemplateAliasMacro(
      vtkImageSquaredDifferenceExecute(
        self, inData0, inData1, stencil,
        inPtr, static_cast<VTK_TT *>(inPtr1), extent, pieceId, output));
    default:
      vtkErrorWithObjectMacro(self, "Execute: Unknown input ScalarType");
    }
}

} // end anonymous namespace

//----------------------------------------------------------------------------
int vtkImageSquaredDifference::RequestData(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // create the thread-local object
  vtkImageSquaredDifferenceTLS tlocal;
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
void vtkImageSquaredDifference::PieceRequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector),
  const int extent[6], vtkIdType pieceId)
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

  int *ext = const_cast<int *>(extent);
  void *inPtr0 = inData0->GetScalarPointerForExtent(ext);
  void *inPtr1 = inData1->GetScalarPointerForExtent(ext);

  vtkImageStencilData *stencil = this->GetStencil();

  switch (inData0->GetScalarType())
    {
    vtkTemplateAliasMacro(
      vtkImageSquaredDifferenceExecute1(
        this, inData0, inData1, stencil,
        static_cast<VTK_TT *>(inPtr0), inPtr1, extent, pieceId,
        &this->ThreadData->Local(pieceId)));
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
    }
}

//----------------------------------------------------------------------------
void vtkImageSquaredDifference::ReduceRequestData(
  vtkInformation *, vtkInformationVector **, vtkInformationVector *)
{
  double sqsum = 0;
  vtkIdType count = 0;

  for (vtkImageSquaredDifferenceTLS::iterator
       iter = this->ThreadData->begin();
       iter != this->ThreadData->end(); ++iter)
    {
    sqsum += iter->SumSquares;
    count += iter->Count;
    }

  if (count == 0)
    {
    count = 1;
    }

  // output values
  this->SquaredDifference = sqsum/count;

  this->SetMinimizable(this->SquaredDifference);
}
