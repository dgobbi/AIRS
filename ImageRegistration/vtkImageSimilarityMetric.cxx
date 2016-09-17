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
  this->Minimizable = 0.0;

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
  os << indent << "Minimizable: " << this->Minimizable << "\n";
}

//----------------------------------------------------------------------------
void vtkImageSimilarityMetric::SetStencilData(vtkImageStencilData *stencil)
{
#if VTK_MAJOR_VERSION >= 6
  this->SetInputData(2, stencil);
#else
  this->SetInput(2, stencil);
#endif
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

  int extent[6] = { 0, -1, 0, -1, 0, -1 };

  bool foundConnection = false;
  int numPorts = ts->Algorithm->GetNumberOfInputPorts();
  for (int inPort = 0; inPort < numPorts; ++inPort)
    {
    int numConnections = ts->Algorithm->GetNumberOfInputConnections(inPort);
    if (numConnections)
      {
      vtkInformation *inInfo =
        ts->InputsInfo[inPort]->GetInformationObject(0);
      vtkImageData *inData = vtkImageData::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));
      if (inData)
        {
        inData->GetExtent(extent);
        foundConnection = true;
        break;
        }
      }
    }

  if (foundConnection)
    {
    // execute the actual method with appropriate extent
    // first find out how many pieces extent can be split into.
    int splitExt[6];
    int total = ts->Algorithm->SplitExtent(
      splitExt, extent, ti->ThreadID, ti->NumberOfThreads);

    if (ti->ThreadID < total &&
        splitExt[1] >= splitExt[0] &&
        splitExt[3] >= splitExt[2] &&
        splitExt[5] >= splitExt[4])
      {
      ts->Algorithm->PieceRequestData(
        ts->Request, ts->InputsInfo, ts->OutputsInfo,
        splitExt, ti->ThreadID);
      }
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
    vtkImageSimilarityMetricThreadStruct *pipelineInfo,
    const int extent[6],
    vtkIdType pieces)
    : PipelineInfo(pipelineInfo),
      NumberOfPieces(pieces)
  {
    for (int i = 0; i < 6; i++)
      {
      this->Extent[i] = extent[i];
      }
  }

  void Initialize() {}
  void operator()(vtkIdType begin, vtkIdType end);
  void Reduce();

private:
  vtkImageSimilarityMetricThreadStruct *PipelineInfo;
  int Extent[6];
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
      splitExt, this->Extent, piece, this->NumberOfPieces);

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
  // start of code copied from vtkThreadedImageAlgorithm

  // setup the threads structure
  vtkImageSimilarityMetricThreadStruct ts;
  ts.Algorithm = this;
  ts.Request = request;
  ts.InputsInfo = inputVector;
  ts.OutputsInfo = outputVector;

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
#if VTK_MAJOR_VERSION >= 6
        this->AllocateOutputData(outData, info, updateExtent);
#else
        this->AllocateOutputData(outData, updateExtent);
#endif
        }
      }
    }

  // copy arrays from first input to output
  int numberOfInputs = this->GetNumberOfInputPorts();
  if (numberOfInputs > 0)
    {
    vtkInformationVector* portInfo = inputVector[0];
    int numberOfConnections = portInfo->GetNumberOfInformationObjects();
    if (numberOfConnections && numberOfOutputs)
      {
      vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
      vtkImageData *inData = vtkImageData::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));
      vtkInformation* outInfo = outputVector->GetInformationObject(0);
      vtkImageData *outData = vtkImageData::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));
      this->CopyAttributeData(inData, outData, inputVector);
      }
    }

  // end of code copied from vtkThreadedImageAlgorithm

#ifdef USE_SMP_THREADED_IMAGE_ALGORITHM
  if (this->EnableSMP)
    {
    // code for vtkSMPTools
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkImageData *inData = vtkImageData::SafeDownCast(
      inInfo->Get(vtkDataObject::DATA_OBJECT()));
    int extent[6];
    inData->GetExtent(extent);

    // do a dummy execution of SplitExtent to compute the number of pieces
    vtkIdType pieces = this->SplitExtent(0, extent, 0, this->NumberOfThreads);

    // create the thread-local object and the functor
    vtkImageSimilarityMetricFunctor functor(&ts, extent, pieces);

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
