/*=========================================================================

  Module: vtkMultiLevelImageRegistration.cxx

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkMultiLevelImageRegistration.h"

#include "vtkImageRegistration.h"

// VTK header files
#include <vtkTimerLog.h>
#include <vtkObjectFactory.h>
#include <vtkImageData.h>
#include <vtkImageStencilData.h>
#include <vtkCommand.h>
#include <vtkMath.h>
#include <vtkDoubleArray.h>
#include <vtkTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkImageSincInterpolator.h>
#include <vtkImageResize.h>
#include <vtkSmartPointer.h>
#include <vtkVersion.h>

// A macro to assist VTK 5 backwards compatibility
#if VTK_MAJOR_VERSION >= 6
#define SET_INPUT_DATA SetInputData
#else
#define SET_INPUT_DATA SetInput
#endif

vtkStandardNewMacro(vtkMultiLevelImageRegistration);

//----------------------------------------------------------------------------
vtkMultiLevelImageRegistration::vtkMultiLevelImageRegistration()
{
  this->Helper = vtkImageRegistration::New();
  this->Level = 0;
  this->NumberOfLevels = 4;
  this->LevelDone = false;

  this->CostValues->Delete();
  this->CostValues = this->Helper->GetCostValues();
  this->CostValues->Register(this);

  this->MetricValues->Delete();
  this->MetricValues = this->Helper->GetMetricValues();
  this->MetricValues->Register(this);

  this->ParameterValues->Delete();
  this->ParameterValues = this->Helper->GetParameterValues();
  this->ParameterValues->Register(this);
}

//----------------------------------------------------------------------------
vtkMultiLevelImageRegistration::~vtkMultiLevelImageRegistration()
{
  if (this->Helper)
    {
    this->Helper->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkMultiLevelImageRegistration::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "Level: " << this->Level << "\n";
  os << indent << "NumberOfLevels: " << this->NumberOfLevels << "\n";
  os << indent << "LevelDone: " << this->LevelDone << "\n";
}

//--------------------------------------------------------------------------
void vtkMultiLevelImageRegistration::Initialize(vtkMatrix4x4 *matrix)
{
  // update our inputs
  this->Update();

  int transformDim = this->TransformDimensionality;
  if (transformDim < 2) { transformDim = 2; }
  if (transformDim > 3) { transformDim = 3; }

  vtkImageData *targetImage = this->GetTargetImage();
  vtkImageData *sourceImage = this->GetSourceImage();

  if (targetImage == NULL || sourceImage == NULL)
    {
    vtkErrorMacro("Initialize: Input images are not set");
    return;
    }

  this->Transform->Identity();
  this->Transform->Concatenate(matrix);

  this->Modified();

  this->InitializeLevel(0);
}

//--------------------------------------------------------------------------
void vtkMultiLevelImageRegistration::InitializeLevel(int level)
{
  this->Level = level;
  this->LevelDone = false;

  // the images to be registered
  vtkImageData *targetImage = this->GetTargetImage();
  vtkImageData *sourceImage = this->GetSourceImage();

  // set up the registration
  vtkImageRegistration *registration = this->Helper;
  registration->CopySettings(this);
  if (level != 0)
    {
    registration->SetInitializerTypeToNone();
    }

  if (level + 1 == this->NumberOfLevels)
    {
    registration->SetSourceImage(sourceImage);
    registration->SetTargetImage(targetImage);
    }
  else
    {
    double blurFactor = 1.0;
    for (int i = level; i < this->NumberOfLevels-1; i++)
      {
      blurFactor *= 2.0;
      }

    // get the pixel spacing for the images
    double targetSpacing[3], sourceSpacing[3];
    targetImage->GetSpacing(targetSpacing);
    sourceImage->GetSpacing(sourceSpacing);

    for (int jj = 0; jj < 3; jj++)
      {
      targetSpacing[jj] = fabs(targetSpacing[jj]);
      sourceSpacing[jj] = fabs(sourceSpacing[jj]);
      }

    double minSpacing = sourceSpacing[0];
    if (minSpacing > sourceSpacing[1])
      {
      minSpacing = sourceSpacing[1];
      }
    if (minSpacing > sourceSpacing[2])
      {
      minSpacing = sourceSpacing[2];
      }

    // reduced resolution: set the blurring
    double spacing[3];
    for (int j = 0; j < 3; j++)
      {
      spacing[j] = blurFactor*minSpacing;
      if (spacing[j] < sourceSpacing[j])
        {
        spacing[j] = sourceSpacing[j];
        }
      }

    // blur source image with Blackman-windowed sinc
    vtkSmartPointer<vtkImageSincInterpolator> sourceBlurKernel =
      vtkSmartPointer<vtkImageSincInterpolator>::New();
    sourceBlurKernel->SetWindowFunctionToBlackman();
    sourceBlurKernel->AntialiasingOn();
    sourceBlurKernel->SetBlurFactors(
      spacing[0]/sourceSpacing[0],
      spacing[1]/sourceSpacing[1],
      spacing[2]/sourceSpacing[2]);

    // reduce the source resolution
    vtkSmartPointer<vtkImageResize> sourceBlur =
      vtkSmartPointer<vtkImageResize>::New();
    sourceBlur->SET_INPUT_DATA(sourceImage);
    sourceBlur->SetResizeMethodToOutputSpacing();
    sourceBlur->SetOutputSpacing(spacing);
    sourceBlur->SetInterpolator(sourceBlurKernel);
    sourceBlur->SetInterpolate((this->InterpolatorType != Nearest));
    sourceBlur->Update();

    // blur target with Blackman-windowed sinc
    vtkSmartPointer<vtkImageSincInterpolator> targetBlurKernel =
      vtkSmartPointer<vtkImageSincInterpolator>::New();
    targetBlurKernel->SetWindowFunctionToBlackman();
    targetBlurKernel->AntialiasingOn();
    targetBlurKernel->SetBlurFactors(
      blurFactor*minSpacing/targetSpacing[0],
      blurFactor*minSpacing/targetSpacing[1],
      blurFactor*minSpacing/targetSpacing[2]);

    // keep target (the image we interpolate) at full resolution
    vtkSmartPointer<vtkImageResize> targetBlur =
      vtkSmartPointer<vtkImageResize>::New();
    targetBlur->SET_INPUT_DATA(targetImage);
    targetBlur->SetResizeMethodToOutputSpacing();
    targetBlur->SetOutputSpacing(targetSpacing);
    targetBlur->SetInterpolator(targetBlurKernel);
    targetBlur->SetInterpolate((this->InterpolatorType != Nearest));
    targetBlur->Update();

    // use blurred images as input to image registration
    registration->SetTransformTolerance(this->TransformTolerance*blurFactor);
    registration->SetSourceImage(sourceBlur->GetOutput());
    registration->SetTargetImage(targetBlur->GetOutput());
    }

  registration->Initialize(this->Transform->GetMatrix());
}

//--------------------------------------------------------------------------
int vtkMultiLevelImageRegistration::ExecuteRegistration()
{
  // reset Abort flag
  this->AbortExecute = 0;
  this->Progress = 0.0;

  this->InvokeEvent(vtkCommand::StartEvent,NULL);
  int converged = 1;

  for (int level = 0; level < this->NumberOfLevels; level++)
    {
    this->Level = level;
    if (level > 0)
      {
      this->InitializeLevel(level);
      }

    this->UpdateProgress(level*1.0/this->NumberOfLevels);

    converged = this->Helper->UpdateRegistration();

    if (this->AbortExecute)
      {
      return 0;
      }

    this->LevelDone = true;
    }

  if (converged && !this->AbortExecute)
    {
    this->UpdateProgress(1.0);
    }

  this->ExecuteTime.Modified();
  this->InvokeEvent(vtkCommand::EndEvent,NULL);

  return converged;
}

//--------------------------------------------------------------------------
int vtkMultiLevelImageRegistration::Iterate()
{
  if (this->LevelDone)
    {
    if (this->Level+1 < this->NumberOfLevels)
      {
      this->InitializeLevel(this->Level + 1);
      }
    else
      {
      return 0;
      }
    }

  vtkImageRegistration *registration = this->Helper;
  int result = registration->Iterate();
  this->Transform->SetMatrix(registration->GetTransform()->GetMatrix());
  this->CostValue = registration->GetCostValue();
  this->NumberOfEvaluations = registration->GetNumberOfEvaluations();

  if (result == 1)
    {
    return 1;
    }

  this->LevelDone = true;

  if (this->Level+1 < this->NumberOfLevels)
    {
    return 1;
    }

  return 0;
}
