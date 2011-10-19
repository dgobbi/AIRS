/*=========================================================================

  Program:   Atamai Classes for VTK
  Module:    $RCSfile: vtkImageRegistration.cxx,v $

Copyright (c) 2005 Atamai, Inc.
All rights reserved.

Use, modification and redistribution of the software, in source or
binary forms, are permitted provided that the following terms and
conditions are met:

1) Redistribution of the source code, in verbatim or modified
   form, must retain the above copyright notice, this license,
   the following disclaimer, and any notices that refer to this
   license and/or the following disclaimer.

2) Redistribution in binary form must include the above copyright
   notice, a copy of this license and the following disclaimer
   in the documentation or with other materials provided with the
   distribution.

3) Modified copies of the source code must be clearly marked as such,
   and must not be misrepresented as verbatim copies of the source code.

THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE SOFTWARE "AS IS"
WITHOUT EXPRESSED OR IMPLIED WARRANTY INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE.  IN NO EVENT SHALL ANY COPYRIGHT HOLDER OR OTHER PARTY WHO MAY
MODIFY AND/OR REDISTRIBUTE THE SOFTWARE UNDER THE TERMS OF THIS LICENSE
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, LOSS OF DATA OR DATA BECOMING INACCURATE
OR LOSS OF PROFIT OR BUSINESS INTERRUPTION) ARISING IN ANY WAY OUT OF
THE USE OR INABILITY TO USE THE SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES.

=========================================================================*/

#include "vtkImageRegistration.h"

// VTK header files
#include "vtkTimerLog.h"
#include "vtkObjectFactory.h"
#include "vtkImageData.h"
#include "vtkImageStencilData.h"
#include "vtkMath.h"
#include "vtkTransform.h"
#include "vtkMatrixToLinearTransform.h"
#include "vtkMatrix4x4.h"
#include "vtkImageReslice.h"
#include "vtkImageShiftScale.h"
#include "vtkCommand.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkAmoebaMinimizer.h"

// Image metric header files
#include "vtkImageMutualInformation.h"
#include "vtkImageCrossCorrelation.h"

// A helper class for the optimizer
struct vtkImageRegistrationInfo
{
  vtkLinearTransform *Transform;
  vtkObject *Optimizer;
  vtkAlgorithm *Metric;

  int TransformType;
  int OptimizerType;
  int MetricType;

  double Center[3];

  int NumberOfEvaluations;
};

//----------------------------------------------------------------------------
vtkImageRegistration* vtkImageRegistration::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret =
    vtkObjectFactory::CreateInstance("vtkImageRegistration");

  if (ret)
    {
    return (vtkImageRegistration*)ret;
    }

  // If the factory was unable to create the object, then create it here.
  return new vtkImageRegistration;
}

//----------------------------------------------------------------------------
vtkImageRegistration::vtkImageRegistration()
{
  this->OptimizerType = vtkImageRegistration::Amoeba;
  this->MetricType = vtkImageRegistration::NormalizedMutualInformation;
  this->InterpolatorType = vtkImageRegistration::Linear;
  this->TransformType = vtkImageRegistration::Rigid;

  this->Transform = vtkTransform::New();
  this->Metric = NULL;
  this->Optimizer = NULL;
  this->Interpolator = NULL;

  this->RegistrationInfo = new vtkImageRegistrationInfo;
  this->RegistrationInfo->Transform = NULL;
  this->RegistrationInfo->Optimizer = NULL;
  this->RegistrationInfo->Metric = NULL;
  this->RegistrationInfo->TransformType = 0;
  this->RegistrationInfo->OptimizerType = 0;
  this->RegistrationInfo->MetricType = 0;
  this->RegistrationInfo->NumberOfEvaluations = 0;

  this->JointHistogramSize[0] = 64;
  this->JointHistogramSize[1] = 64;

  this->ImageReslice = vtkImageReslice::New();
  this->TargetImageQuantizer = vtkImageShiftScale::New();
  this->SourceImageQuantizer = vtkImageShiftScale::New();

  this->Value = 0.0;

  this->MetricTolerance = 1e-4;
  this->TransformTolerance = 1e-1;
  this->MaximumNumberOfIterations = 500;

  // we have the image inputs and the optional stencil input
  this->SetNumberOfInputPorts(3);
  this->SetNumberOfOutputPorts(0);
}

//----------------------------------------------------------------------------
vtkImageRegistration::~vtkImageRegistration()
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

  if (this->RegistrationInfo)
    {
    delete this->RegistrationInfo;
    }

  if (this->ImageReslice)
    {
    this->ImageReslice->Delete();
    }
  if (this->SourceImageQuantizer)
    {
    this->SourceImageQuantizer->Delete();
    }
  if (this->TargetImageQuantizer)
    {
    this->TargetImageQuantizer->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkImageRegistration::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "OptimizerType: " << this->OptimizerType << "\n";
  os << indent << "MetricType: " << this->MetricType << "\n";
  os << indent << "InterpolatorType: " << this->InterpolatorType << "\n";
  os << indent << "TransformType: " << this->TransformType << "\n";
  os << indent << "MetricTolerance: " << this->MetricTolerance << "\n";
  os << indent << "TransformTolerance: " << this->TransformTolerance << "\n";
  os << indent << "MaximumNumberOfIterations: "
     << this->MaximumNumberOfIterations << "\n";
  os << indent << "JointHistogramSize: " << this->JointHistogramSize[0] << " "
     << this->JointHistogramSize[1] << "\n";
  os << indent << "Value: " << this->Value << "\n";
  os << indent << "NumberOfEvaluations: "
     << this->RegistrationInfo->NumberOfEvaluations << "\n";
}

//----------------------------------------------------------------------------
int vtkImageRegistration::GetNumberOfEvaluations()
{
  return this->RegistrationInfo->NumberOfEvaluations;
}

//----------------------------------------------------------------------------
void vtkImageRegistration::SetTargetImage(vtkImageData *input)
{
  // Ask the superclass to connect the input.
  this->SetNthInputConnection(0, 0, (input ? input->GetProducerPort() : 0));
}

//----------------------------------------------------------------------------
vtkImageData* vtkImageRegistration::GetTargetImage()
{
  if (this->GetNumberOfInputConnections(0) < 1)
    {
    return NULL;
    }
  return vtkImageData::SafeDownCast(this->GetExecutive()->GetInputData(0, 0));
}

//----------------------------------------------------------------------------
void vtkImageRegistration::SetSourceImage(vtkImageData *input)
{
  // Ask the superclass to connect the input.
  this->SetNthInputConnection(0, 1, (input ? input->GetProducerPort() : 0));
}

//----------------------------------------------------------------------------
vtkImageData* vtkImageRegistration::GetSourceImage()
{
  if (this->GetNumberOfInputConnections(1) < 1)
    {
    return NULL;
    }
  return vtkImageData::SafeDownCast(this->GetExecutive()->GetInputData(1, 0));
}

//----------------------------------------------------------------------------
void vtkImageRegistration::SetTargetImageStencil(vtkImageStencilData *stencil)
{
  // if stencil is null, then set the input port to null
  this->SetNthInputConnection(2, 0,
    (stencil ? stencil->GetProducerPort() : 0));
}

//----------------------------------------------------------------------------
vtkImageStencilData* vtkImageRegistration::GetTargetImageStencil()
{
  if (this->GetNumberOfInputConnections(2) < 1)
    {
    return NULL;
    }
  return vtkImageStencilData::SafeDownCast(
    this->GetExecutive()->GetInputData(2, 0));
}

//--------------------------------------------------------------------------
namespace {

void vtkSetTransformParameters(vtkImageRegistrationInfo *registrationInfo)
{
  vtkAmoebaMinimizer* optimizer =
    vtkAmoebaMinimizer::SafeDownCast(registrationInfo->Optimizer);
  vtkTransform* transform =
    vtkTransform::SafeDownCast(registrationInfo->Transform);

  double tx = optimizer->GetParameterValue(0);
  double ty = optimizer->GetParameterValue(1);
  double tz = optimizer->GetParameterValue(2);

  double rx = optimizer->GetParameterValue(3);
  double ry = optimizer->GetParameterValue(4);
  double rz = optimizer->GetParameterValue(5);

  double s = 1.0;
  if (registrationInfo->TransformType != vtkImageRegistration::Rigid)
    {
    s = optimizer->GetParameterValue(6);
    }

  double *center = registrationInfo->Center;

  // compute quaternion from rotation parameters
  double qs = rx*rx + ry*ry + rz*rz;
  double qc = 1.0 - qs;
  while (qc < 0)
    {
    // for rotations past 180 degrees
    qs -= 1.0;
    qc = 1.0 - qs;
    rx = -rx;
    ry = -ry;
    rz = -rz;
    }
  qs = sqrt(qs);
  qc = sqrt(qc);
  double theta = atan2(qs,qc)*180/vtkMath::DoublePi();
  if (qs > 0)
    {
    rx /= qs;
    ry /= qs;
    rz /= qs;
    }

  transform->Identity();
  transform->PostMultiply();
  transform->Translate(-center[0], -center[1], -center[2]);
  transform->RotateWXYZ(theta, rx, ry, rz);
  transform->Scale(s,s,s);
  transform->Translate(center[0], center[1], center[2]);
  transform->Translate(tx,ty,tz);
}

//--------------------------------------------------------------------------
void vtkEvaluateFunction(void * arg)
{
  vtkImageRegistrationInfo *registrationInfo =
    static_cast<vtkImageRegistrationInfo*>(arg);

  double val = 0.0;

  vtkAmoebaMinimizer* optimizer =
    vtkAmoebaMinimizer::SafeDownCast(registrationInfo->Optimizer);
  vtkTransform* transform =
    vtkTransform::SafeDownCast(registrationInfo->Transform);
  vtkImageMutualInformation *miMetric =
    vtkImageMutualInformation::SafeDownCast(registrationInfo->Metric);
  vtkImageCrossCorrelation *ccMetric =
    vtkImageCrossCorrelation::SafeDownCast(registrationInfo->Metric);

  vtkSetTransformParameters(registrationInfo);

  registrationInfo->Metric->Update();

  switch (registrationInfo->MetricType)
    {
    case vtkImageRegistration::CrossCorrelation:
      val = - ccMetric->GetCrossCorrelation();
      break;
    case vtkImageRegistration::NormalizedCrossCorrelation:
      val = - ccMetric->GetNormalizedCrossCorrelation();
      break;
    case vtkImageRegistration::MutualInformation:
      val = - miMetric->GetMutualInformation();
      break;
    case vtkImageRegistration::NormalizedMutualInformation:
      val = - miMetric->GetNormalizedMutualInformation();
      break;
    }

  optimizer->SetFunctionValue(val);

  registrationInfo->NumberOfEvaluations++;
}

} // end anonymous namespace

//--------------------------------------------------------------------------
void vtkImageRegistration::ComputeImageRange(
  vtkImageData *data, double range[2])
{
  data->GetScalarRange(range);

  //range[0] = 0.0;
  if (range[0] >= range[1])
    {
    range[1] = range[0] + 1.0;
    }
}

//--------------------------------------------------------------------------
void vtkImageRegistration::Initialize(vtkMatrix4x4 *matrix)
{
  // update our inputs
  this->Update();

  vtkImageData *targetImage = this->GetTargetImage();
  vtkImageData *sourceImage = this->GetSourceImage();

  if (targetImage == NULL || sourceImage == NULL)
    {
    vtkErrorMacro("Initialize: Input images are not set");
    return;
    }

  // get the target image center
  double bounds[6];
  double center[3];
  double size[3];
  targetImage->GetBounds(bounds);
  center[0] = 0.5*(bounds[0] + bounds[1]);
  center[1] = 0.5*(bounds[2] + bounds[3]);
  center[2] = 0.5*(bounds[4] + bounds[5]);
  size[0] = (bounds[1] - bounds[0]);
  size[1] = (bounds[3] - bounds[2]);
  size[2] = (bounds[5] - bounds[4]);

  vtkTransform *transform =
    vtkTransform::SafeDownCast(this->Transform);

  transform->Identity();

  if (matrix)
    {
    // if a matrix was given, use it to make a baseline
    // for the registration transform
    vtkTransform *initialTransform = vtkTransform::New();
    initialTransform->Concatenate(matrix);
    transform->SetInput(initialTransform);
    initialTransform->Delete();
    }
  else
    {
    // if no matrix was given, set an initial translation
    // from one image center to the other image center
    double sbounds[6];
    double scenter[3];
    sourceImage->GetBounds(sbounds);
    scenter[0] = 0.5*(sbounds[0] + sbounds[1]);
    scenter[1] = 0.5*(sbounds[2] + sbounds[3]);
    scenter[2] = 0.5*(sbounds[4] + sbounds[5]);
    vtkTransform *initialTransform = vtkTransform::New();
    initialTransform->Translate(
      scenter[0]-center[0], scenter[1]-center[1], scenter[2]-center[2]);
    transform->SetInput(initialTransform);
    initialTransform->Delete();
    }

  if ((this->MetricType ==
       vtkImageRegistration::MutualInformation ||
       this->MetricType ==
       vtkImageRegistration::NormalizedMutualInformation) &&
      this->InterpolatorType == vtkImageRegistration::Nearest &&
      this->JointHistogramSize[0] <= 256 &&
      this->JointHistogramSize[1] <= 256)
    {
    // If nearest-neighbor interpolation is used, then the image instensity
    // can be quantized during initialization, instead of being done at each
    // iteration of the registration.

    double sourceImageRange[2];
    this->ComputeImageRange(sourceImage, sourceImageRange);

    double sourceScale = ((this->JointHistogramSize[1] - 1)/
      (sourceImageRange[1] - sourceImageRange[0]));
    // The "0.5/sourceScale" causes the vtkImageShiftScale filter to
    // round the value, instead of truncating it.
    double sourceShift = (-sourceImageRange[0] + 0.5/sourceScale);

    vtkImageShiftScale *sourceQuantizer = this->SourceImageQuantizer;
    sourceQuantizer->SetInput(sourceImage);
    sourceQuantizer->SetOutputScalarTypeToUnsignedChar();
    sourceQuantizer->ClampOverflowOn();
    sourceQuantizer->SetShift(sourceShift);
    sourceQuantizer->SetScale(sourceScale);
    sourceQuantizer->Update();
    sourceImage = sourceQuantizer->GetOutput();

    double targetImageRange[2];
    this->ComputeImageRange(targetImage, targetImageRange);

    double targetScale = ((this->JointHistogramSize[0] - 1)/
      (targetImageRange[1] - targetImageRange[0]));
    double targetShift = (-targetImageRange[0] + 0.5/targetScale);

    vtkImageShiftScale *targetQuantizer = this->TargetImageQuantizer;
    targetQuantizer->SetInput(targetImage);
    targetQuantizer->SetOutputScalarTypeToUnsignedChar();
    targetQuantizer->ClampOverflowOn();
    targetQuantizer->SetShift(targetShift);
    targetQuantizer->SetScale(targetScale);
    targetQuantizer->Update();
    targetImage = targetQuantizer->GetOutput();
    }

  vtkImageReslice *reslice = this->ImageReslice;
  reslice->SetInput(sourceImage);
  reslice->SetInformationInput(targetImage);
  reslice->SetStencil(this->GetTargetImageStencil());
  reslice->SetResliceTransform(this->Transform);
  reslice->GenerateStencilOutputOn();
  switch (this->InterpolatorType)
    {
    case vtkImageRegistration::Nearest:
      reslice->SetInterpolationModeToNearestNeighbor();
      break;
    case vtkImageRegistration::Linear:
      reslice->SetInterpolationModeToLinear();
      break;
    case vtkImageRegistration::Cubic:
      reslice->SetInterpolationModeToCubic();
      break;
    }

  if (this->Metric)
    {
    this->Metric->RemoveAllInputs();
    this->Metric->Delete();
    }

  switch (this->MetricType)
    {
    case vtkImageRegistration::CrossCorrelation:
    case vtkImageRegistration::NormalizedCrossCorrelation:
      {
      vtkImageCrossCorrelation *metric = vtkImageCrossCorrelation::New();
      this->Metric = metric;

      metric->SetInput(targetImage);
      metric->SetInputConnection(1, reslice->GetOutputPort());
      metric->SetInputConnection(2, reslice->GetStencilOutputPort());
      }
      break;

    case vtkImageRegistration::MutualInformation:
    case vtkImageRegistration::NormalizedMutualInformation:
      {
      vtkImageMutualInformation *metric = vtkImageMutualInformation::New();
      this->Metric = metric;

      metric->SetInput(targetImage);
      metric->SetInputConnection(1, reslice->GetOutputPort());
      metric->SetInputConnection(2, reslice->GetStencilOutputPort());
      metric->SetNumberOfBins(this->JointHistogramSize);

      double targetImageRange[2], sourceImageRange[2];
      this->ComputeImageRange(targetImage, targetImageRange);
      this->ComputeImageRange(sourceImage, sourceImageRange);

      metric->SetBinOrigin(targetImageRange[0], sourceImageRange[0]);
      metric->SetBinSpacing(
        (targetImageRange[1] - targetImageRange[0])/
          (this->JointHistogramSize[0]-1),
        (sourceImageRange[1] - sourceImageRange[0])/
          (this->JointHistogramSize[1]-1));
      }
      break;
    }

  if (this->Optimizer != NULL)
    {
    this->Optimizer->Delete();
    }

  vtkAmoebaMinimizer *optimizer = vtkAmoebaMinimizer::New();
  this->Optimizer = optimizer;
  optimizer->SetTolerance(this->MetricTolerance);
  optimizer->SetParameterTolerance(this->TransformTolerance);
  optimizer->SetMaxIterations(this->MaximumNumberOfIterations);

  this->RegistrationInfo->Transform = this->Transform;
  this->RegistrationInfo->Optimizer = this->Optimizer;
  this->RegistrationInfo->Metric = this->Metric;

  this->RegistrationInfo->TransformType = this->TransformType;
  this->RegistrationInfo->OptimizerType = this->OptimizerType;
  this->RegistrationInfo->MetricType = this->MetricType;

  this->RegistrationInfo->NumberOfEvaluations = 0;

  this->RegistrationInfo->Center[0] = center[0];
  this->RegistrationInfo->Center[1] = center[1];
  this->RegistrationInfo->Center[2] = center[2];

  // use golden ratio for amoeba
  optimizer->SetExpansionRatio(1.618);
  optimizer->SetContractionRatio(0.618);
  optimizer->SetFunction(&vtkEvaluateFunction,
                         (void*)(this->RegistrationInfo));

  /*
  double rotation[3], translation[3], scale3[3];
  transform->PreMultiply();
  transform->Translate(center[0],center[1],center[2]);
  transform->PostMultiply();
  transform->Translate(-center[0],-center[1],-center[2]);
  transform->GetOrientation(rotation);
  transform->GetScale(scale3);
  transform->GetPosition(translation);
  double scale = pow(fabs(scale3[0]*scale3[1]*scale3[2]), 1.0/3);
  */

  // compute minimum spacing of target image
  double spacing[3];
  targetImage->GetSpacing(spacing);
  double minspacing = spacing[0];
  minspacing = ((minspacing < spacing[1]) ? minspacing : spacing[1]);
  minspacing = ((minspacing < spacing[2]) ? minspacing : spacing[2]);

  // compute a radius of gyration
  double r = 1.0;
  int dims = 0;
  for (int ii = 0; ii < 3; ii++)
    {
    if (size[ii] > 1e-6)
      {
      r *= size[ii];
      dims++;
      }
    }
  r = 0.41*pow(r, 1.0/dims);

  // compute parameter scales
  double tscale = this->TransformTolerance*10;
  tscale = ((tscale >= minspacing) ? tscale : minspacing);
  double rscale = tscale/r;
  double sscale = tscale/r;
  if (rscale > 0.5)
    {
    rscale = 0.5;
    }
  if (sscale > 0.1)
    {
    sscale = 0.1;
    }

  optimizer->Initialize();

  optimizer->SetParameterValue(0, 0);
  optimizer->SetParameterScale(0, tscale);
  optimizer->SetParameterValue(1, 0);
  optimizer->SetParameterScale(1, tscale);
  optimizer->SetParameterValue(2, 0);
  optimizer->SetParameterScale(2, tscale);

  optimizer->SetParameterValue(3, 0);
  optimizer->SetParameterScale(3, rscale);
  optimizer->SetParameterValue(4, 0);
  optimizer->SetParameterScale(4, rscale);
  optimizer->SetParameterValue(5, 0);
  optimizer->SetParameterScale(5, rscale);

  if (this->TransformType != vtkImageRegistration::Rigid)
    {
    optimizer->SetParameterValue(6, 1);
    optimizer->SetParameterScale(6, sscale);
    }

  this->Modified();
}

//--------------------------------------------------------------------------
int vtkImageRegistration::ExecuteRegistration()
{
  // reset Abort flag
  this->AbortExecute = 0;
  this->Progress = 0.0;

  this->InvokeEvent(vtkCommand::StartEvent,NULL);

  int converged = 0;

  vtkAmoebaMinimizer *optimizer =
    vtkAmoebaMinimizer::SafeDownCast(this->Optimizer);

  if (optimizer)
    {
    int n = this->MaximumNumberOfIterations;
    if (n <= 0)
      {
      n = VTK_INT_MAX;
      }
    for (int i = 0; i < n && !converged; i++)
      {
      this->UpdateProgress(i*1.0/n);
      if (this->AbortExecute)
        {
        break;
        }
      converged = !optimizer->Iterate();
      vtkSetTransformParameters(this->RegistrationInfo);
      this->Value = optimizer->GetFunctionValue();
      }

    if (converged && !this->AbortExecute)
      {
      this->UpdateProgress(1.0);
      }
    }

  this->ExecuteTime.Modified();
  this->InvokeEvent(vtkCommand::EndEvent,NULL);

  return converged;
}

//--------------------------------------------------------------------------
int vtkImageRegistration::Iterate()
{
  vtkAmoebaMinimizer *optimizer =
    vtkAmoebaMinimizer::SafeDownCast(this->Optimizer);

  if (optimizer)
    {
    int result = optimizer->Iterate();
    if (optimizer->GetIterations() >= this->MaximumNumberOfIterations)
      {
      result = 0;
      }
    vtkSetTransformParameters(this->RegistrationInfo);
    this->Value = optimizer->GetFunctionValue();
    return result;
    }

  return 0;
}

//--------------------------------------------------------------------------
int vtkImageRegistration::UpdateRegistration()
{
  this->Update();
  return this->ExecuteRegistration();
}

//----------------------------------------------------------------------------
int vtkImageRegistration::FillInputPortInformation(int port,
                                                   vtkInformation* info)
{
  if (port == 2)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageStencilData");
    // the stencil input is optional
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }
  else
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    }

  return 1;
}

//----------------------------------------------------------------------------
int vtkImageRegistration::FillOutputPortInformation(
  int port, vtkInformation* info)
{
  return 1;
}

//----------------------------------------------------------------------------
int vtkImageRegistration::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector))
{
  return 1;
}

//----------------------------------------------------------------------------
int vtkImageRegistration::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector))
{
  int inExt[6];

  // source image
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), inExt);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inExt, 6);

  // target image
  inInfo = inputVector[1]->GetInformationObject(0);
  inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), inExt);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inExt, 6);

  // stencil for target image
  if (this->GetNumberOfInputConnections(2) > 0)
    {
    vtkInformation *inInfo2 = inputVector[2]->GetInformationObject(0);
    inInfo2->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inExt, 6);
    }

  return 1;
}
//----------------------------------------------------------------------------
int vtkImageRegistration::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector))
{
  return 1;
}

//----------------------------------------------------------------------------
int vtkImageRegistration::ProcessRequest(vtkInformation* request,
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
