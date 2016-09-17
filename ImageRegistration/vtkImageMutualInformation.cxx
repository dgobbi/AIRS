/*=========================================================================

  Module: vtkImageMutualInformation.cxx

  Copyright (c) 2006 Atamai, Inc.
  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageMutualInformation.h"

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

vtkStandardNewMacro(vtkImageMutualInformation);

//----------------------------------------------------------------------------
// Data needed for each thread.
class vtkImageMutualInformationThreadData
{
public:
  vtkImageMutualInformationThreadData() : Data(0) {}

  vtkIdType *Data;
};

class vtkImageMutualInformationTLS
  : public vtkImageSimilarityMetricTLS<vtkImageMutualInformationThreadData>
{
};

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageMutualInformation::vtkImageMutualInformation()
{
  this->NumberOfBins[0] = 64;
  this->NumberOfBins[1] = 64;

  this->BinOrigin[0] = 0.0;
  this->BinOrigin[1] = 0.0;

  this->BinSpacing[0] = 1.0;
  this->BinSpacing[1] = 1.0;

  this->OutputScalarType = VTK_FLOAT;

  this->MutualInformation = 0.0;
  this->NormalizedMutualInformation = 0.0;

  this->ThreadData = 0;

  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkImageMutualInformation::~vtkImageMutualInformation()
{
}

//----------------------------------------------------------------------------
void vtkImageMutualInformation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "NumberOfBins: " << this->NumberOfBins[0] << " "
     << this->NumberOfBins[1] << "\n";
  os << indent << "BinOrigin: " << this->BinOrigin[0] << " "
     << this->BinOrigin[1] << "\n";
  os << indent << "BinSpacing: " << this->BinSpacing[0] << " "
     << this->BinSpacing[1] << "\n";

  os << indent << "MutualInformation: " << this->MutualInformation << "\n";
  os << indent << "NormalizedMutualInformation: "
     << this->NormalizedMutualInformation << "\n";
}

//----------------------------------------------------------------------------
int vtkImageMutualInformation::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  int outWholeExt[6];
  double outOrigin[3];
  double outSpacing[3];

  outWholeExt[0] = 0;
  outWholeExt[1] = this->NumberOfBins[0] - 1;
  outWholeExt[2] = 0;
  outWholeExt[3] = this->NumberOfBins[1] - 1;
  outWholeExt[4] = 0;
  outWholeExt[5] = 0;

  outOrigin[0] = this->BinOrigin[0];
  outOrigin[1] = this->BinOrigin[1];
  outOrigin[2] = 0.0;

  outSpacing[0] = this->BinSpacing[0];
  outSpacing[1] = this->BinSpacing[1];
  outSpacing[2] = 1.0;

  int outScalarType = this->OutputScalarType;
  int outComponents = 1;

  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),outWholeExt,6);
  outInfo->Set(vtkDataObject::ORIGIN(), outOrigin, 3);
  outInfo->Set(vtkDataObject::SPACING(), outSpacing, 3);

  vtkDataObject::SetPointDataActiveScalarInfo(
      outInfo, outScalarType, outComponents);

  return 1;
}

//----------------------------------------------------------------------------
// anonymous namespace for internal functions
namespace {

//----------------------------------------------------------------------------
template<class T1, class T2>
void vtkImageMutualInformationExecute(
  vtkImageMutualInformation *self,
  vtkImageData *inData0, vtkImageData *inData1, vtkImageStencilData *stencil,
  T1 *inPtr, T2 *inPtr1, int extent[6],
  vtkIdType *outPtr, int numBins[2], double binOrigin[2], double binSpacing[2],
  vtkIdType pieceId)
{
  vtkImageStencilIterator<T1>
    inIter(inData0, stencil, extent, ((pieceId == 0) ? self : NULL));
  vtkImageStencilIterator<T2>
    inIter1(inData1, stencil, extent, NULL);

  int pixelInc = inData0->GetNumberOfScalarComponents();
  int pixelInc1 = inData1->GetNumberOfScalarComponents();

  static double xmin = 0;
  static double ymin = 0;
  double xmax = numBins[0] - 1;
  double ymax = numBins[1] - 1;
  double xshift = -binOrigin[0];
  double yshift = -binOrigin[1];
  double xscale = 1.0/binSpacing[0];
  double yscale = 1.0/binSpacing[1];
  int outIncY = numBins[1];

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

        x += xshift;
        x *= xscale;

        y += yshift;
        y *= yscale;

        x = (x > xmin ? x : xmin);
        x = (x < xmax ? x : xmax);

        y = (y > ymin ? y : ymin);
        y = (y < ymax ? y : ymax);

        int xi = static_cast<int>(x + 0.5);
        int yi = static_cast<int>(y + 0.5);

        vtkIdType *outPtr1 = outPtr + yi*outIncY + xi;

        (*outPtr1)++;

        inPtr += pixelInc;
        inPtr1 += pixelInc1;
        }
      }
    inIter.NextSpan();
    inIter1.NextSpan();
    }
}

//----------------------------------------------------------------------------
void vtkImageMutualInformationExecutePreScaled(
  vtkImageMutualInformation *self,
  vtkImageData *inData0, vtkImageData *inData1, vtkImageStencilData *stencil,
  unsigned char *inPtr, unsigned char *inPtr1, int extent[6],
  vtkIdType *outPtr, int numBins[2], vtkIdType pieceId)
{
  vtkImageStencilIterator<unsigned char>
    inIter(inData0, stencil, extent, ((pieceId == 0) ? self : NULL));

  vtkImageStencilIterator<unsigned char>
    inIter1(inData1, stencil, extent, NULL);

  int pixelInc = inData0->GetNumberOfScalarComponents();
  int pixelInc1 = inData1->GetNumberOfScalarComponents();

  int xmax = numBins[0] - 1;
  int ymax = numBins[1] - 1;
  int outIncY = numBins[1];

  // iterate over all spans in the stencil
  while (!inIter.IsAtEnd())
    {
    if (inIter.IsInStencil())
      {
      inPtr = inIter.BeginSpan();
      unsigned char *inPtrEnd = inIter.EndSpan();
      inPtr1 = inIter1.BeginSpan();

      // iterate over all voxels in the span
      while (inPtr != inPtrEnd)
        {
        int x = *inPtr;
        int y = *inPtr1;

        x = (x < xmax ? x : xmax);
        y = (y < ymax ? y : ymax);

        vtkIdType *outPtr1 = outPtr + y*outIncY + x;

        (*outPtr1)++;

        inPtr += pixelInc;
        inPtr1 += pixelInc1;
        }
      }
    inIter.NextSpan();
    inIter1.NextSpan();
    }
}

//----------------------------------------------------------------------------
// copy one row of the joint histogram to the output, with conversion
// but without type range checking
template<class T>
void vtkImageMutualInformationCopyRow(
  vtkIdType *xyHist, T *outPtr, int outStart, int outEnd)
{
  int n = outEnd - outStart + 1;
  xyHist += outStart;

  do
    {
    *outPtr++ = static_cast<T>(*xyHist++);
    }
  while (--n);
}

} // end anonymous namespace

//----------------------------------------------------------------------------
int vtkImageMutualInformation::RequestData(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // create the thread-local object
  vtkImageMutualInformationTLS tlocal;
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
template<class T1>
void vtkImageMutualInformationExecute1(
  vtkImageMutualInformation *self,
  vtkImageData *inData0, vtkImageData *inData1, vtkImageStencilData *stencil,
  T1 *inPtr, void *inPtr1, int extent[6],
  vtkIdType *outPtr, int numBins[2], double binOrigin[2], double binSpacing[2],
  vtkIdType pieceId)
{
  switch (inData1->GetScalarType())
    {
    vtkTemplateAliasMacro(
      vtkImageMutualInformationExecute(
        self, inData0, inData1, stencil,
        inPtr, static_cast<VTK_TT *>(inPtr1), extent,
        outPtr, numBins, binOrigin, binSpacing, pieceId));
    default:
      vtkErrorWithObjectMacro(self, "Execute: Unknown input ScalarType");
    }
}

} // end anonymous namespace

//----------------------------------------------------------------------------
// This method is passed a input and output region, and executes the filter
// algorithm to fill the output from the input.
// It just executes a switch statement to call the correct function for
// the regions data types.
void vtkImageMutualInformation::PieceRequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector),
  const int pieceExtent[6], vtkIdType pieceId)
{
  vtkImageMutualInformationThreadData *threadLocal =
    &this->ThreadData->Local(pieceId);

  vtkIdType *outPtr = threadLocal->Data;

  if (outPtr == 0)
    {
    // initialize the joint histogram to zero
    vtkIdType outIncY = this->NumberOfBins[0];
    vtkIdType outCount = outIncY;
    outCount *= this->NumberOfBins[1];

    threadLocal->Data = new vtkIdType[outCount];
    outPtr = threadLocal->Data;

    vtkIdType *outPtr1 = outPtr;
    do { *outPtr1++ = 0; } while (--outCount > 0);
    }

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

  int extent[6];
  for (int i = 0; i < 6; i += 2)
    {
    int j = i + 1;
    extent[i] = pieceExtent[i];
    extent[i] = ((extent[i] > inExt0[i]) ? extent[i] : inExt0[i]);
    extent[i] = ((extent[i] > inExt1[i]) ? extent[i] : inExt1[i]);
    extent[j] = pieceExtent[j];
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

  double *binOrigin = this->BinOrigin;
  double *binSpacing = this->BinSpacing;
  int *numBins = this->NumberOfBins;
  int maxX = numBins[0] - 1;
  int maxY = numBins[1] - 1;

  if (vtkMath::Floor(binOrigin[0] + 0.5) == 0 &&
      vtkMath::Floor(binOrigin[1] + 0.5) == 0 &&
      vtkMath::Floor(binOrigin[0] + binSpacing[0]*maxX + 0.5) == maxX &&
      vtkMath::Floor(binOrigin[1] + binSpacing[1]*maxY + 0.5) == maxY &&
      inData0->GetScalarType() == VTK_UNSIGNED_CHAR &&
      inData1->GetScalarType() == VTK_UNSIGNED_CHAR)
    {
    vtkImageMutualInformationExecutePreScaled(
      this, inData0, inData1, stencil,
      static_cast<unsigned char *>(inPtr0),
      static_cast<unsigned char *>(inPtr1),
      extent, outPtr, numBins, pieceId);
    }
  else switch (inData0->GetScalarType())
    {
    vtkTemplateAliasMacro(
      vtkImageMutualInformationExecute1(
        this, inData0, inData1, stencil,
        static_cast<VTK_TT *>(inPtr0), inPtr1,
        extent, outPtr, numBins, binOrigin, binSpacing,
        pieceId));
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
    }
}

//----------------------------------------------------------------------------
void vtkImageMutualInformation::ReduceRequestData(
  vtkInformation *, vtkInformationVector **,
  vtkInformationVector *outputVector)
{
  // get the output
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  int updateExtent[6];
  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
               updateExtent);
  vtkImageData *outData = vtkImageData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // get the output pointer
  void *outPtr = outData->GetScalarPointerForExtent(updateExtent);
  int outScalarType = outData->GetScalarType();
  int outScalarSize = outData->GetScalarSize();

  // get the dimensions of the joint histogram
  int nx = this->NumberOfBins[0];
  int ny = this->NumberOfBins[1];

  // various variables for computing the mutual information
  double xEntropy = 0;
  double yEntropy = 0;
  double xyEntropy = 0;

  // allocate space to accumulate results
  vtkIdType *xyHist = new vtkIdType[2*nx + ny];
  vtkIdType *xHist = xyHist + nx;

  // clear xHist to zero
  int ix;
  for (ix = 0; ix < nx; ++ix)
    {
    xHist[ix] = 0;
    }

  // piece together the joint histogram results from each thread
  for (int iy = 0; iy < ny; ++iy)
    {
    int outStart = updateExtent[0];
    int outEnd = updateExtent[1];
    if (iy < updateExtent[2] || iy > updateExtent[3])
      {
      outStart = 0;
      outEnd = -1;
      }

    // clear the output
    for (ix = 0; ix < nx; ++ix)
      {
      xyHist[ix] = 0;
      }

    // add the contribution from thread j
    vtkIdType a = 0;
    for (vtkImageMutualInformationTLS::iterator
         iter = this->ThreadData->begin();
         iter != this->ThreadData->end(); ++iter)
      {
      if (iter->Data)
        {
        vtkIdType *outPtr2 = iter->Data + static_cast<vtkIdType>(nx)*iy;

        for (ix = 0; ix < nx; ++ix)
          {
          vtkIdType c = *outPtr2++;
          xyHist[ix] += c;
          a += c;
          }
        }
      }

    // copy this row of the joint histogram to the output
    if (outStart <= outEnd)
      {
      switch (outScalarType)
        {
        vtkTemplateAliasMacro(
          vtkImageMutualInformationCopyRow(
            xyHist, static_cast<VTK_TT *>(outPtr), outStart, outEnd));
        default:
          vtkErrorMacro("Execute: Unknown output ScalarType");
        }
      // increment outPtr to the next row
      outPtr =
        static_cast<char *>(outPtr) + outScalarSize*(outEnd - outStart + 1);
      }

    // compute the entropy of second image
    double da = static_cast<double>(a);
    if (da > 0)
      {
      yEntropy += da*log(da);
      }

    // compute joint entropy
    for (ix = 0; ix < nx; ++ix)
      {
      vtkIdType c = xyHist[ix];
      xHist[ix] += c;
      double dc = static_cast<double>(c);
      if (dc > 0)
        {
        xyEntropy += dc*log(dc);
        }
      }
    }

  // compute total pixel count and entropy of first image
  vtkIdType count = 0;
  for (ix = 0; ix < nx; ++ix)
    {
    vtkIdType b = xHist[ix];
    count += b;
    double db = static_cast<double>(b);
    if (db > 0)
      {
      xEntropy += db*log(db);
      }
    }

  for (vtkImageMutualInformationTLS::iterator
       iter = this->ThreadData->begin();
       iter != this->ThreadData->end(); ++iter)
    {
    // delete the temporary memory
    delete [] iter->Data;
    }

  delete [] xyHist;

  // minimum possible values
  double mutualInformation = 0.0;
  double normalizedMutualInformation = 1.0;

  if (count)
    {
    // correct for total voxel count, convert to negative
    double dc = static_cast<double>(count);
    double ldc = log(dc);
    xEntropy = -xEntropy/dc + ldc;
    yEntropy = -yEntropy/dc + ldc;
    xyEntropy = -xyEntropy/dc + ldc;

    // compute the mutual information
    mutualInformation = xEntropy + yEntropy - xyEntropy;

    // compute the normalized mutual information (Studholme 1999)
    normalizedMutualInformation = (xEntropy + yEntropy)/xyEntropy;
    }

  this->MutualInformation = mutualInformation;
  this->NormalizedMutualInformation = normalizedMutualInformation;

  this->SetMinimizable(mutualInformation);
}
