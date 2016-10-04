/*=========================================================================

  Module: vtkImageRegistrationBase.cxx

  Copyright (c) 2006 Atamai, Inc.
  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkImageRegistrationBase.h"

// VTK header files
#include <vtkImageData.h>
#include <vtkImageStencilData.h>
#include <vtkMath.h>
#include <vtkDoubleArray.h>
#include <vtkTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkCommand.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkImageHistogramStatistics.h>
#include <vtkVersion.h>

// Interpolator header files
#include "vtkLabelInterpolator.h"

// Optimizer header files
#include "vtkNelderMeadMinimizer.h"
#include "vtkPowellMinimizer.h"

// Image metric header files
#include "vtkImageSquaredDifference.h"
#include "vtkImageMutualInformation.h"
#include "vtkImageCorrelationRatio.h"
#include "vtkImageCrossCorrelation.h"
#include "vtkImageNeighborhoodCorrelation.h"

// C header files
#include <math.h>

// A macro to assist VTK 5 backwards compatibility
#if VTK_MAJOR_VERSION >= 6
#define SET_INPUT_DATA SetInputData
#define SET_STENCIL_DATA SetStencilData
#else
#define SET_INPUT_DATA SetInput
#define SET_STENCIL_DATA SetStencil
#endif

//----------------------------------------------------------------------------
vtkImageRegistrationBase::vtkImageRegistrationBase()
{
  this->OptimizerType = vtkImageRegistrationBase::Powell;
  this->MetricType = vtkImageRegistrationBase::MutualInformation;
  this->InterpolatorType = vtkImageRegistrationBase::Linear;
  this->TransformType = vtkImageRegistrationBase::Rigid;
  this->InitializerType = vtkImageRegistrationBase::None;
  this->TransformDimensionality = 3;

  this->Transform = vtkTransform::New();
  this->Metric = NULL;
  this->Optimizer = NULL;
  this->Interpolator = NULL;

  this->JointHistogramSize[0] = 64;
  this->JointHistogramSize[1] = 64;
  this->SourceImageRange[0] = 0.0;
  this->SourceImageRange[1] = -1.0;
  this->TargetImageRange[0] = 0.0;
  this->TargetImageRange[1] = -1.0;

  this->CostValue = 0.0;
  this->NumberOfEvaluations = 0;

  this->CollectValues = false;
  this->MetricValues = vtkDoubleArray::New();
  this->CostValues = vtkDoubleArray::New();
  this->ParameterValues = vtkDoubleArray::New();

  this->CostTolerance = 1e-4;
  this->TransformTolerance = 1e-1;
  this->MaximumNumberOfIterations = 500;
  this->MaximumNumberOfEvaluations = 5000;

  // we have the image inputs and the optional stencil input
  this->SetNumberOfInputPorts(4);
  this->SetNumberOfOutputPorts(0);
}

//----------------------------------------------------------------------------
vtkImageRegistrationBase::~vtkImageRegistrationBase()
{
  // delete vtk objects
  if (this->Optimizer)
    {
    this->Optimizer->Delete();
    }
  if (this->Metric)
    {
    this->Metric->Delete();
    }
  if (this->Interpolator)
    {
    this->Interpolator->Delete();
    }
  if (this->Transform)
    {
    this->Transform->Delete();
    }
  if (this->MetricValues)
    {
    this->MetricValues->Delete();
    }
  if (this->CostValues)
    {
    this->CostValues->Delete();
    }
  if (this->ParameterValues)
    {
    this->ParameterValues->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkImageRegistrationBase::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "OptimizerType: " << this->OptimizerType << "\n";
  os << indent << "MetricType: " << this->MetricType << "\n";
  os << indent << "InterpolatorType: " << this->InterpolatorType << "\n";
  os << indent << "TransformType: " << this->TransformType << "\n";
  os << indent << "TransformDimensionality: "
     << this->TransformDimensionality << "\n";
  os << indent << "InitializerType: " << this->InitializerType << "\n";
  os << indent << "CostTolerance: " << this->CostTolerance << "\n";
  os << indent << "TransformTolerance: " << this->TransformTolerance << "\n";
  os << indent << "MaximumNumberOfIterations: "
     << this->MaximumNumberOfIterations << "\n";
  os << indent << "MaximumNumberOfEvaluations: "
     << this->MaximumNumberOfEvaluations << "\n";
  os << indent << "JointHistogramSize: " << this->JointHistogramSize[0] << " "
     << this->JointHistogramSize[1] << "\n";
  os << indent << "SourceImageRange: " << this->SourceImageRange[0] << " "
     << this->SourceImageRange[1] << "\n";
  os << indent << "TargetImageRange: " << this->TargetImageRange[0] << " "
     << this->TargetImageRange[1] << "\n";
  os << indent << "CostValue: " << this->CostValue << "\n";
  os << indent << "CollectValues: "
     << (this->CollectValues ? "On\n" : "Off\n");
  os << indent << "MetricValues: " << this->MetricValues << "\n";
  os << indent << "CostValues: " << this->CostValues << "\n";
  os << indent << "ParameterValues: " << this->ParameterValues << "\n";
  os << indent << "NumberOfEvaluations: "
     << this->NumberOfEvaluations << "\n";
}

//----------------------------------------------------------------------------
vtkLinearTransform *vtkImageRegistrationBase::GetTransform()
{
  return this->Transform;
}

//----------------------------------------------------------------------------
int vtkImageRegistrationBase::GetNumberOfEvaluations()
{
  return this->NumberOfEvaluations;
}

//----------------------------------------------------------------------------
void vtkImageRegistrationBase::CopySettings(vtkImageRegistrationBase *other)
{
  bool modified = false;

  if (this->MetricType != other->MetricType)
    {
    this->MetricType = other->MetricType;
    modified = true;
    }
  if (this->OptimizerType != other->OptimizerType)
    {
    this->OptimizerType = other->OptimizerType;
    modified = true;
    }
  if (this->InterpolatorType != other->InterpolatorType)
    {
    this->InterpolatorType = other->InterpolatorType;
    modified = true;
    }
  if (this->TransformType != other->TransformType)
    {
    this->TransformType = other->TransformType;
    modified = true;
    }
  if (this->TransformDimensionality != other->TransformDimensionality)
    {
    this->TransformDimensionality = other->TransformDimensionality;
    modified = true;
    }
  if (this->InitializerType != other->InitializerType)
    {
    this->InitializerType = other->InitializerType;
    modified = true;
    }
  if (this->JointHistogramSize[0] != other->JointHistogramSize[0] ||
      this->JointHistogramSize[1] != other->JointHistogramSize[1])
    {
    this->JointHistogramSize[0] = other->JointHistogramSize[0];
    this->JointHistogramSize[1] = other->JointHistogramSize[1];
    modified = true;
    }
  if (this->SourceImageRange[0] != other->SourceImageRange[0] ||
      this->SourceImageRange[1] != other->SourceImageRange[1])
    {
    this->SourceImageRange[0] = other->SourceImageRange[0];
    this->SourceImageRange[1] = other->SourceImageRange[1];
    modified = true;
    }
  if (this->TargetImageRange[0] != other->TargetImageRange[0] ||
      this->TargetImageRange[1] != other->TargetImageRange[1])
    {
    this->TargetImageRange[0] = other->TargetImageRange[0];
    this->TargetImageRange[1] = other->TargetImageRange[1];
    modified = true;
    }
  if (this->CostTolerance != other->CostTolerance)
    {
    this->CostTolerance = other->CostTolerance;
    modified = true;
    }
  if (this->TransformTolerance != other->TransformTolerance)
    {
    this->TransformTolerance = other->TransformTolerance;
    modified = true;
    }
  if (this->MaximumNumberOfIterations != other->MaximumNumberOfIterations)
    {
    this->MaximumNumberOfIterations = other->MaximumNumberOfIterations;
    modified = true;
    }
  if (this->MaximumNumberOfEvaluations != other->MaximumNumberOfEvaluations)
    {
    this->MaximumNumberOfEvaluations = other->MaximumNumberOfEvaluations;
    modified = true;
    }
  if (this->CollectValues != other->CollectValues)
    {
    this->CollectValues = other->CollectValues;
    modified = true;
    }
}

//----------------------------------------------------------------------------
void vtkImageRegistrationBase::SetTargetImage(vtkImageData *input)
{
  // Ask the superclass to connect the input.
#if VTK_MAJOR_VERSION >= 6
  this->SetInputDataInternal(1, input);
#else
  this->SetNthInputConnection(1, 0, (input ? input->GetProducerPort() : 0));
#endif
}

//----------------------------------------------------------------------------
vtkImageData* vtkImageRegistrationBase::GetTargetImage()
{
  if (this->GetNumberOfInputConnections(0) < 1)
    {
    return NULL;
    }
  return vtkImageData::SafeDownCast(this->GetExecutive()->GetInputData(1, 0));
}

//----------------------------------------------------------------------------
void vtkImageRegistrationBase::SetSourceImage(vtkImageData *input)
{
  // Ask the superclass to connect the input.
#if VTK_MAJOR_VERSION >= 6
  this->SetInputDataInternal(0, input);
#else
  this->SetNthInputConnection(0, 0, (input ? input->GetProducerPort() : 0));
#endif
}

//----------------------------------------------------------------------------
vtkImageData* vtkImageRegistrationBase::GetSourceImage()
{
  if (this->GetNumberOfInputConnections(1) < 1)
    {
    return NULL;
    }
  return vtkImageData::SafeDownCast(this->GetExecutive()->GetInputData(0, 0));
}

//----------------------------------------------------------------------------
void vtkImageRegistrationBase::SetSourceImageStencil(vtkImageStencilData *stencil)
{
  // if stencil is null, then set the input port to null
#if VTK_MAJOR_VERSION >= 6
  this->SetInputDataInternal(2, stencil);
#else
  this->SetNthInputConnection(2, 0,
    (stencil ? stencil->GetProducerPort() : 0));
#endif
}

//----------------------------------------------------------------------------
vtkImageStencilData* vtkImageRegistrationBase::GetSourceImageStencil()
{
  if (this->GetNumberOfInputConnections(2) < 1)
    {
    return NULL;
    }
  return vtkImageStencilData::SafeDownCast(
    this->GetExecutive()->GetInputData(2, 0));
}

//----------------------------------------------------------------------------
void vtkImageRegistrationBase::SetTargetImageStencil(
  vtkImageStencilData *stencil)
{
  // if stencil is null, then set the input port to null
#if VTK_MAJOR_VERSION >= 6
  this->SetInputDataInternal(3, stencil);
#else
  this->SetNthInputConnection(3, 0,
    (stencil ? stencil->GetProducerPort() : 0));
#endif
}

//----------------------------------------------------------------------------
vtkImageStencilData* vtkImageRegistrationBase::GetTargetImageStencil()
{
  if (this->GetNumberOfInputConnections(3) < 1)
    {
    return NULL;
    }
  return vtkImageStencilData::SafeDownCast(
    this->GetExecutive()->GetInputData(3, 0));
}

//--------------------------------------------------------------------------
void vtkImageRegistrationBase::ComputeImageRange(
  vtkImageData *data, vtkImageStencilData *stencil, double range[2])
{
  vtkImageHistogramStatistics *hist =
    vtkImageHistogramStatistics::New();
  hist->SET_STENCIL_DATA(stencil);
  hist->SET_INPUT_DATA(data);
  hist->SetActiveComponent(0);
  hist->Update();

  range[0] = hist->GetMinimum();
  range[1] = hist->GetMaximum();

  if (range[0] >= range[1])
    {
    range[1] = range[0] + 1.0;
    }

  hist->SET_INPUT_DATA(NULL);
  hist->Delete();
}

//--------------------------------------------------------------------------
int vtkImageRegistrationBase::UpdateRegistration()
{
  this->Update();
  return this->ExecuteRegistration();
}

//----------------------------------------------------------------------------
int vtkImageRegistrationBase::FillInputPortInformation(int port,
                                                   vtkInformation* info)
{
  if (port == 2 || port == 3)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageStencilData");
    // the stencil input is optional
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }
  else
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    }

  return 1;
}

//----------------------------------------------------------------------------
int vtkImageRegistrationBase::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* vtkNotUsed(info))
{
  return 1;
}

//----------------------------------------------------------------------------
int vtkImageRegistrationBase::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *vtkNotUsed(outputVector))
{
  return 1;
}

//----------------------------------------------------------------------------
int vtkImageRegistrationBase::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector))
{
  int inExt[6];

  // source image
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), inExt);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inExt, 6);

  // stencil for source image
  if (this->GetNumberOfInputConnections(2) > 0)
    {
    vtkInformation *inInfo2 = inputVector[2]->GetInformationObject(0);
    inInfo2->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inExt, 6);
    }

  // target image
  inInfo = inputVector[1]->GetInformationObject(0);
  inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), inExt);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inExt, 6);

  // stencil for target image
  if (this->GetNumberOfInputConnections(3) > 0)
    {
    vtkInformation *inInfo3 = inputVector[3]->GetInformationObject(0);
    inInfo3->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inExt, 6);
    }

  return 1;
}
//----------------------------------------------------------------------------
int vtkImageRegistrationBase::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *vtkNotUsed(outputVector))
{
  return 1;
}

//----------------------------------------------------------------------------
int vtkImageRegistrationBase::ProcessRequest(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // generate the data oject
  if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_OBJECT()))
    {
    return 1;
    }
  // generate the data
  if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
    {
    return this->RequestData(request, inputVector, outputVector);
    }

  // execute information
  if (request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
    {
    return this->RequestInformation(request, inputVector, outputVector);
    }

  // propagate update extent
  if (request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
    {
    return this->RequestUpdateExtent(request, inputVector, outputVector);
    }

  return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}
