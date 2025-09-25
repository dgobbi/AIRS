/*=========================================================================

  Module: vtkImageSimilarityMetric.cxx

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageSimilarityMetric.h"

#include <vtkImageData.h>
#include <vtkImageStencilData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkMultiThreader.h>
#include <vtkVersion.h>

#include "vtkImageSimilarityMetricInternals.h"

#include <math.h>

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageSimilarityMetric::vtkImageSimilarityMetric()
{
  this->InputRange[0][0] = 0.0;
  this->InputRange[0][1] = 0.0;
  this->InputRange[1][0] = 0.0;
  this->InputRange[1][1] = 0.0;

  this->Value = 0.0;
  this->Cost = 0.0;

  this->SetNumberOfInputPorts(3);
  this->SetNumberOfOutputPorts(0);
}

//----------------------------------------------------------------------------
vtkImageSimilarityMetric::~vtkImageSimilarityMetric()
{
}

//----------------------------------------------------------------------------
void vtkImageSimilarityMetric::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Stencil: " << this->GetStencil() << "\n";
  os << indent << "InputRange: ("
     << this->InputRange[0][0] << ", " << this->InputRange[0][1] << "), ("
     << this->InputRange[1][0] << ", " << this->InputRange[1][1] << ")\n";
  os << indent << "Value: " << this->Value << "\n";
  os << indent << "Cost: " << this->Cost << "\n";
}

//----------------------------------------------------------------------------
void vtkImageSimilarityMetric::SetStencilData(vtkImageStencilData *stencil)
{
  this->SetInputData(2, stencil);
}

//----------------------------------------------------------------------------
vtkImageStencilData *vtkImageSimilarityMetric::GetStencil()
{
  if (this->GetNumberOfInputConnections(2) < 1)
  {
    return NULL;
  }
  return vtkImageStencilData::SafeDownCast(
    this->GetExecutive()->GetInputData(2, 0));
}

//----------------------------------------------------------------------------
void vtkImageSimilarityMetric::SetInputRange(int i, const double r[2])
{
  if (i >= 0 && i < 2)
  {
    if (r[0] != this->InputRange[i][0] ||
        r[1] != this->InputRange[i][1])
    {
      this->InputRange[i][0] = r[0];
      this->InputRange[i][1] = r[1];
      this->Modified();
    }
  }
}

//----------------------------------------------------------------------------
void vtkImageSimilarityMetric::GetInputRange(int i, double r[2])
{
  if (i >= 0 && i < 2)
  {
    r[0] = this->InputRange[i][0];
    r[1] = this->InputRange[i][1];
  }
  else
  {
    r[0] = 0.0;
    r[1] = 0.0;
  }
}

//----------------------------------------------------------------------------
int vtkImageSimilarityMetric::FillInputPortInformation(
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
int vtkImageSimilarityMetric::RequestUpdateExtent(
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

//----------------------------------------------------------------------------
struct vtkImageSimilarityMetricThreadStruct
{
  static VTK_THREAD_RETURN_TYPE ThreadExecute(void *arg);

  vtkImageSimilarityMetric *Algorithm;
  vtkInformation *Request;
  vtkInformationVector **InputsInfo;
  vtkInformationVector *OutputsInfo;
  int Extent[6];
};

//----------------------------------------------------------------------------
// override from vtkThreadedImageAlgorithm to split input extent, instead
// of splitting the output extent
VTK_THREAD_RETURN_TYPE
vtkImageSimilarityMetricThreadStruct::ThreadExecute(void *arg)
{
  vtkMultiThreader::ThreadInfo *ti =
    static_cast<vtkMultiThreader::ThreadInfo *>(arg);
  vtkImageSimilarityMetricThreadStruct *ts =
    static_cast<vtkImageSimilarityMetricThreadStruct *>(ti->UserData);

  // execute the actual method with appropriate extent
  // first find out how many pieces extent can be split into.
  int splitExt[6];
  int total = ts->Algorithm->SplitExtent(
    splitExt, ts->Extent, ti->ThreadID, ti->NumberOfThreads);

  if (ti->ThreadID < total &&
      splitExt[1] >= splitExt[0] &&
      splitExt[3] >= splitExt[2] &&
      splitExt[5] >= splitExt[4])
  {
    ts->Algorithm->PieceRequestData(
      ts->Request, ts->InputsInfo, ts->OutputsInfo,
      splitExt, ti->ThreadID);
  }

  return VTK_THREAD_RETURN_VALUE;
}

//----------------------------------------------------------------------------
#ifdef USE_SMP_THREADED_IMAGE_ALGORITHM
// Functor for vtkSMPTools execution
class vtkImageSimilarityMetricFunctor
{
public:
  // Create the functor, provide all info it needs to execute.
  vtkImageSimilarityMetricFunctor(
    vtkImageSimilarityMetricThreadStruct *pipelineInfo, vtkIdType pieces)
    : PipelineInfo(pipelineInfo), NumberOfPieces(pieces) {}

  void Initialize() {}
  void operator()(vtkIdType begin, vtkIdType end);
  void Reduce();

private:
  vtkImageSimilarityMetricThreadStruct *PipelineInfo;
  vtkIdType NumberOfPieces;
};

// Called by vtkSMPTools to execute the algorithm over specific pieces.
void vtkImageSimilarityMetricFunctor::operator()(
  vtkIdType begin, vtkIdType end)
{
  vtkImageSimilarityMetricThreadStruct *ts = this->PipelineInfo;

  for (vtkIdType piece = begin; piece < end; piece++)
  {
    int splitExt[6] = { 0, -1, 0, -1, 0, -1 };

    vtkIdType total = ts->Algorithm->SplitExtent(
      splitExt, ts->Extent, piece, this->NumberOfPieces);

    // check for valid piece and extent
    if (piece < total &&
        splitExt[0] <= splitExt[1] &&
        splitExt[2] <= splitExt[3] &&
        splitExt[4] <= splitExt[5])
    {
      ts->Algorithm->PieceRequestData(
        ts->Request, ts->InputsInfo, ts->OutputsInfo, splitExt, piece);
    }
  }
}

// Called by vtkSMPTools once the multi-threading has finished.
void vtkImageSimilarityMetricFunctor::Reduce()
{
  vtkImageSimilarityMetric *self =
    static_cast<vtkImageSimilarityMetric *>(this->PipelineInfo->Algorithm);
  vtkImageSimilarityMetricThreadStruct *ts = this->PipelineInfo;

  self->ReduceRequestData(ts->Request, ts->InputsInfo, ts->OutputsInfo);
}
#endif

//----------------------------------------------------------------------------
// override from vtkThreadedImageAlgorithm to customize the multithreading
int vtkImageSimilarityMetric::RequestData(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // setup the threads structure
  vtkImageSimilarityMetricThreadStruct ts;
  ts.Algorithm = this;
  ts.Request = request;
  ts.InputsInfo = inputVector;
  ts.OutputsInfo = outputVector;

  vtkInformation *inInfo0 = inputVector[0]->GetInformationObject(0);
  vtkInformation *inInfo1 = inputVector[1]->GetInformationObject(0);

  vtkImageData *inData0 = vtkImageData::SafeDownCast(
    inInfo0->Get(vtkDataObject::DATA_OBJECT()));
  vtkImageData *inData1 = vtkImageData::SafeDownCast(
    inInfo1->Get(vtkDataObject::DATA_OBJECT()));

  // allocate the output data
  int numberOfOutputs = this->GetNumberOfOutputPorts();
  if (numberOfOutputs > 0)
  {
    for (int i = 0; i < numberOfOutputs; ++i)
    {
      vtkInformation* info = outputVector->GetInformationObject(i);
      vtkImageData *outData = vtkImageData::SafeDownCast(
        info->Get(vtkDataObject::DATA_OBJECT()));
      if (outData)
      {
        int updateExtent[6];
        info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
                  updateExtent);
        this->AllocateOutputData(outData, info, updateExtent);
      }
      // copy arrays from first input to output
      if (i == 0)
      {
        this->CopyAttributeData(inData0, outData, inputVector);
      }
    }
  }

  // Get the intersection of the input extents
  int inExt1[6];
  inData0->GetExtent(ts.Extent);
  inData1->GetExtent(inExt1);

  for (int i = 0; i < 6; i += 2)
  {
    int j = i + 1;
    ts.Extent[i] = ((ts.Extent[i] > inExt1[i]) ? ts.Extent[i] : inExt1[i]);
    ts.Extent[j] = ((ts.Extent[j] < inExt1[j]) ? ts.Extent[j] : inExt1[j]);
    if (ts.Extent[i] > ts.Extent[j])
    {
      // no overlap, nothing to do!
      return 1;
    }
  }

#ifdef USE_SMP_THREADED_IMAGE_ALGORITHM
  if (this->EnableSMP)
  {
    // code for vtkSMPTools

    // do a dummy execution of SplitExtent to compute the number of pieces
    vtkIdType pieces = this->SplitExtent(
      0, ts.Extent, 0, this->NumberOfThreads);

    // create the functor
    vtkImageSimilarityMetricFunctor functor(&ts, pieces);

    bool debug = this->Debug;
    this->Debug = false;
    vtkSMPTools::For(0, pieces, functor);
    this->Debug = debug;
  }
  else
#endif
  {
    // code for vtkMultiThreader
    this->Threader->SetNumberOfThreads(this->NumberOfThreads);
    this->Threader->SetSingleMethod(
      vtkImageSimilarityMetricThreadStruct::ThreadExecute, &ts);

    // always shut off debugging to avoid threading problems with GetMacros
    int debug = this->Debug;
    this->Debug = 0;
    this->Threader->SingleMethodExecute();
    this->Debug = debug;

    this->ReduceRequestData(request, inputVector, outputVector);
  }

  return 1;
}
