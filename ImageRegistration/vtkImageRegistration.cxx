/*=========================================================================

  Program:   Atamai Classes for VTK
  Module:    $RCSfile: vtkImageRegistration.cxx,v $
  Creator:   Kevin Wang <kwang@atamai.com>
  Language:  C++
  Author:    $Author: kwang $
  Date:      $Date: 2006/08/28 18:16:20 $
  Version:   $Revision: 1.2 $

==========================================================================

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
#include "vtkMatrix4x4.h"
#include "vtkHomogeneousTransform.h"
#include "vtkMatrixToHomogeneousTransform.h"
#include "vtkImageReslice.h"
#include "vtkCommand.h"
#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION <= 4)
#else
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#endif

// Optimizer header files
#include "vtkAmoebaMinimizer.h"

// Image metric header files
#include "vtkCalcCrossCorrelation.h"
#include "vtkImageMutualInformation.h"

// Transform header files
#include "vtkCenteredTransform.h"

//--------------------------------------------------------------------------
// optimizer names and optimizer parameter names

const char * vtkOptimizerNames[] = {
  "vtkAmoebaMinimizer", 
  0
};

const char * vtkAmoebaMinimizerParameterNames[] = {
  "MaxNumberOfIterations",
  "Tolerance",
  0
};

//--------------------------------------------------------------------------
// metric names and metric parameter names

const char * vtkMetricNames[] = {
  "vtkImageCrossCorrelation",
  "vtkImageMutualInfomation", 
  0
};

const char * vtkCalcCrossCorrelationParameterNames[] = {
  0
};

const char * vtkImageMutualInformationParameterNames[] = {
  "AComponentExtent",
  "BComponentExtent",
  0          
};

//--------------------------------------------------------------------------
// transform names and tranform parameter names

const char * vtkTransformNames[] = {
  "vtkCenteredTransform",
  0
};

const char * vtkCenteredTransformParameterNames[] = {
  "CenterX",
  "CenterY",
  "CenterZ",
  "TranslationX", 
  "TranslationY",
  "TranslationZ",
  "RotationAngleX",
  "RotationAngleY",
  "RotationAngleZ",
  "IsotropicScale",
  "TranslationScale",
  "RotationScale",
  "ScalingScale",
  0
};

//--------------------------------------------------------------------------
// interpolator names

const char * vtkInterpolatorNames[] = {
  "vtkNearestNeighborInterpolator",
  "vtkLinearInterpolator",
  "vtkCubicInterpolator", 
  0
};

//----------------------------------------------------------------------------
vtkImageRegistration* 
vtkImageRegistration::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = 
    vtkObjectFactory::CreateInstance("vtkImageRegistration");

  if(ret)
    {
    return (vtkImageRegistration*)ret;
    }

  // If the factory was unable to create the object, then create it here.
  return new vtkImageRegistration;
}

//----------------------------------------------------------------------------
vtkImageRegistration::vtkImageRegistration()
{
  this->OptimizerType             = VTK_TYPE_INVALID;
  this->MetricType                = VTK_TYPE_INVALID;
  this->InterpolatorType          = VTK_TYPE_INVALID;
  this->TransformType             = VTK_TYPE_INVALID;

  this->Transform                 = NULL; 
  this->Interpolator              = NULL; 
  this->Metric                    = NULL; 
  this->Optimizer                 = NULL; 

  LastTransform                   = vtkMatrixToHomogeneousTransform::New();

#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION < 4)
  this->OptimizerParameters       = std::vector<double>();
  this->MetricParameters          = std::vector<double>();
  this->TransformParameters       = std::vector<double>();
#else
  this->OptimizerParameters       = vtkstd::vector<double>();
  this->MetricParameters          = vtkstd::vector<double>();
  this->TransformParameters       = vtkstd::vector<double>();
#endif

  this->Value                     = 0.0;
  this->CurrentIteration          = 0;
  this->CurPosition[0]            = 0.0;
  this->CurPosition[1]            = 0.0;
  this->CurPosition[2]            = 0.0;
  this->CurPosition[3]            = 0.0;
  this->CurPosition[4]            = 0.0;
  this->CurPosition[5]            = 0.0;
  this->CurPosition[6]            = 0.0;
  this->CurPosition[7]            = 0.0;
  this->CurPosition[8]            = 0.0;
  this->CurPosition[9]            = 0.0;
  this->CurPosition[10]           = 0.0;
  this->CurPosition[11]           = 0.0;
}

//----------------------------------------------------------------------------
vtkImageRegistration::~vtkImageRegistration() 
{
  // delete vtk objects
  if ( this->LastTransform != NULL )
    {
    LastTransform->Delete();
    }
}
  
//----------------------------------------------------------------------------
void 
vtkImageRegistration::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkProcessObject::PrintSelf(os,indent);
  
  os << indent << "Optimizer: " << this->OptimizerType << "\n";
  os << indent << "Metric: " << this->MetricType << "\n";
  os << indent << "Interpolator: " << this->InterpolatorType << "\n";
  os << indent << "Transform: " << this->TransformType << "\n";
}  

//----------------------------------------------------------------------------
void 
vtkImageRegistration::SetFixedImage(vtkImageData *input)
{
#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION <= 4)
  this->vtkProcessObject::SetNthInput(0, input);
#else
  if (input)
    {
    this->SetNthInputConnection(0, 0, input->GetProducerPort());
    }
  else
    {
    this->SetNthInputConnection(0, 0, 0);
    }
#endif
}

//----------------------------------------------------------------------------
vtkImageData*
vtkImageRegistration::GetFixedImage()
{
#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION <= 4)
  if (this->GetNumberOfInputs() < 1)    
    {
    return NULL;
    }
  return (vtkImageData *)(this->Inputs[0]);
#else
  if (this->GetNumberOfInputConnections(0) < 1)
    {
    return NULL;
    }
  return vtkImageData::SafeDownCast(this->GetExecutive()->GetInputData(0, 0));
#endif
}

//----------------------------------------------------------------------------
void 
vtkImageRegistration::SetMovingImage(vtkImageData *input)
{
#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION <= 4)
  this->vtkProcessObject::SetNthInput(1, input);
#else
  if (input)
    {
    this->SetNthInputConnection(0, 1, input->GetProducerPort());
    }
  else
    {
    this->SetNthInputConnection(0, 1, 0);
    }
#endif
}

//----------------------------------------------------------------------------
vtkImageData*
vtkImageRegistration::GetMovingImage()
{
#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION <= 4)
  if (this->GetNumberOfInputs() < 2)    
    {
    return NULL;
    }
  return (vtkImageData *)(this->Inputs[1]);
#else
  if (this->GetNumberOfInputConnections(0) < 2)
    {
    return NULL;
    }
  return vtkImageData::SafeDownCast(this->GetExecutive()->GetInputData(0, 1));
#endif
}

//----------------------------------------------------------------------------
void 
vtkImageRegistration::SetFixedImageStencil(vtkImageStencilData *stencil)
{
#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION <= 4)
  this->vtkProcessObject::SetNthInput(2, stencil);
#else
  if (stencil)
    {
    this->SetNthInputConnection(0, 2, stencil->GetProducerPort());
    }
  else
    {
    this->SetNthInputConnection(0, 2, 0);
    }
#endif
}

//----------------------------------------------------------------------------
vtkImageStencilData*
vtkImageRegistration::GetFixedImageStencil()
{
#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION <= 4)
  if (this->NumberOfInputs < 3)
    {
    return NULL;
    }
  else
    {
    return (vtkImageStencilData *)(this->Inputs[2]);
    }
#else
  if (this->GetNumberOfInputConnections(0) < 3)
    {
    return NULL;
    }
  return vtkImageStencilData::SafeDownCast(this->GetExecutive()->GetInputData(0, 2));
#endif
}

//----------------------------------------------------------------------------
const char*
vtkImageRegistration::GetOptimizerName(int i)
{
  if ( i < 0 || i > VTK_NUMBER_OF_OPTIMIZERS )
    {
    return NULL;
    }
 
  return vtkOptimizerNames[i];
}

//----------------------------------------------------------------------------
const char* 
vtkImageRegistration::GetOptimizerParameterName(int optimizer, int parameter)
{
  int i, numberOfParameters;
  const char ** tempNames = NULL;

  if ( optimizer < 0 || optimizer > VTK_NUMBER_OF_OPTIMIZERS )
    {
    return NULL;
    }
  
  switch ( optimizer )
    {
    case VTK_OPTIMIZER_AMOEBA:
      tempNames = vtkAmoebaMinimizerParameterNames;
      break;
    }     

  numberOfParameters = 0;
  for ( i = 0; tempNames[i] != 0; i++)
    {
    numberOfParameters++;
    }
  
  if ( parameter < 0 || parameter > numberOfParameters )
    {
    return NULL;
    }
    
  return tempNames[parameter];
}

//----------------------------------------------------------------------------
void 
vtkImageRegistration::SetOptimizerParameter(const char * name, double value)
{
  unsigned int i, numberOfParameters;
  const char ** tempNames = NULL;

  if ( this->OptimizerType < 0 ||  this->OptimizerType 
       >= VTK_NUMBER_OF_OPTIMIZERS )
    {
    vtkErrorMacro( << "SetOptmizerParameter: set optimizer type before set optimizer parameters");
    return ;
    }
  
  switch ( this->OptimizerType )
    {
    case VTK_OPTIMIZER_AMOEBA:
      tempNames = vtkAmoebaMinimizerParameterNames;
      break;
    }     

  numberOfParameters = 0;
  for ( i = 0; tempNames[i] != 0; i++ )
    {
    numberOfParameters++;
    }
  
  if ( this->OptimizerParameters.size() != numberOfParameters ) 
    {
    this->OptimizerParameters.resize( numberOfParameters );
    }
  
  for ( i = 0; i < numberOfParameters; i++ )
    {
    if ( strcmp( name, tempNames[i] ) == 0 )
      {
      if ( this->OptimizerParameters[i] != value )
        {
        this->OptimizerParameters[i] = value;
        this->Modified();
        }
      }
    } 
}

//----------------------------------------------------------------------------
const char*
vtkImageRegistration::GetMetricName(int i)
{
  if ( i < 0 || i > VTK_NUMBER_OF_METRICS )
    {
    return NULL;
    }

  return vtkMetricNames[i];
}

//----------------------------------------------------------------------------
const char*
vtkImageRegistration::GetMetricParameterName(int metric, int parameter)
{
  int i, numberOfParameters;
  const char ** tempNames = NULL;

  if ( metric < 0 || metric > VTK_NUMBER_OF_METRICS )
    {
    return NULL;
    }
  
  switch ( metric )
    {
    case VTK_METRIC_NORMALIZED_CROSS_CORRELATION:
      tempNames = vtkCalcCrossCorrelationParameterNames;
      break;
    case VTK_METRIC_NORMALIZED_MUTUAL_INFORMATION:
      tempNames = vtkImageMutualInformationParameterNames;
      break;
    }     

  numberOfParameters = 0;
  for ( i = 0; tempNames[i] != 0; i++)
    {
    numberOfParameters++;
    }
  
  if ( parameter < 0 || parameter > numberOfParameters )
    {
    return NULL;
    }
    
  return tempNames[parameter];
}

//----------------------------------------------------------------------------
void 
vtkImageRegistration::SetMetricParameter(const char * name , double value)
{
  int i, numberOfParameters;
  const char ** tempNames = NULL;

  if ( this->MetricType < 0 ||
       this->MetricType >= VTK_NUMBER_OF_METRICS )
    {
    vtkErrorMacro(<< "SetMetricParameter: set metric type before set metric parameters");
    return;
    }

  switch ( this->MetricType )
    {
    case VTK_METRIC_NORMALIZED_CROSS_CORRELATION:
      tempNames = vtkCalcCrossCorrelationParameterNames;
      break;
    case VTK_METRIC_NORMALIZED_MUTUAL_INFORMATION:
      tempNames = vtkImageMutualInformationParameterNames;
      break;
    default:
      vtkErrorMacro(<< "SetMetricParameter: set metric type before set metric parameters");
    }     

  numberOfParameters = 0;
  for ( i = 0; tempNames[i] != 0; i++)
    {
    numberOfParameters++;
    }
  
  if ( this->MetricParameters.size() !=
       static_cast<unsigned int>(numberOfParameters) )
    {
    this->MetricParameters.resize( numberOfParameters );
    }
  
  for ( i = 0; i < numberOfParameters; i++ )
    {
    if ( strcmp( name, tempNames[i] ) == 0 )
      {
      if ( this->MetricParameters[i] != value )
        {
        this->MetricParameters[i] = value;
        this->Modified();
        }
      }
    }
}

//----------------------------------------------------------------------------
const char*
vtkImageRegistration::GetInterpolatorName(int i)
{
  if ( i < 0 || i > VTK_NUMBER_OF_INTERPOLATORS )
    {
    return NULL;
    }

  return vtkInterpolatorNames[i];
}

//----------------------------------------------------------------------------
const char*
vtkImageRegistration::GetTransformName(int i)
{
  if ( i < 0 || i > VTK_NUMBER_OF_TRANSFORMS )
    {
    return NULL;
    }

  return vtkTransformNames[i];
}

const char *
vtkImageRegistration::GetTransformParameterName(int transform, int parameter)
{
  int i, numberOfParameters;
  const char ** tempNames = NULL;

  if ( transform < 0 || transform > VTK_NUMBER_OF_TRANSFORMS )
    {
    return NULL;
    }

  switch ( transform )
    {
    case VTK_TRANSFORM_CENTERED:
      tempNames = vtkCenteredTransformParameterNames;
      break;
    }     

  numberOfParameters = 0;
  for ( i = 0; tempNames[i] != 0; i++)
    {
    numberOfParameters++;
    }
  
  if ( parameter < 0 || parameter > numberOfParameters )
    {
    return NULL;
    }
    
  return tempNames[parameter];
}


//----------------------------------------------------------------------------
void 
vtkImageRegistration::SetTransformParameter(const char * name, double value)
{
  int i, numberOfParameters;
  const char ** tempNames = NULL;

  if ( this->TransformType < 0 ||  this->TransformType 
       >= VTK_NUMBER_OF_TRANSFORMS )
    {
    vtkErrorMacro(<< "SetTransformParameter: set transform type before set transform parameters");
    return;
    }

  switch ( this->TransformType )
    {
    case VTK_TRANSFORM_CENTERED:
      tempNames = vtkCenteredTransformParameterNames;
      break;
    default:
      vtkErrorMacro(<< "SetTransformParameter: Unknown TransformType");
    }     

  numberOfParameters = 0;
  for ( i = 0; tempNames[i] != 0; i++)
    {
    numberOfParameters++;
    }
  
  if ( this->TransformParameters.size() !=
       static_cast<unsigned int>(numberOfParameters) )
    {
    this->TransformParameters.resize( numberOfParameters );
    }

  for ( i = 0; i < numberOfParameters; i++ )
    {
    if ( strcmp( name, tempNames[i] ) == 0 )
      {
      if ( this->TransformParameters[i] != value )
        {
        this->TransformParameters[i] = value;
        this->Modified();
        }
      }
    }
}

//----------------------------------------------------------------------------
double 
vtkImageRegistration::GetTransformParameter(const char * name)
{
  int i, numberOfParameters;
  double value = 0.0;
  const char ** tempNames = NULL;

  if ( this->TransformType < 0 ||
       this->TransformType >= VTK_NUMBER_OF_TRANSFORMS )
    {
    return 0.0;
    }

  switch ( this->TransformType )
    {
    case VTK_TRANSFORM_CENTERED:
      tempNames = vtkCenteredTransformParameterNames;
      break;
    default:
      vtkErrorMacro(<< "GetLastTransformParameter: Unknown TransformType");
    }     

  numberOfParameters = 0;
  for ( i = 0; tempNames[i] != 0; i++)
    {
    numberOfParameters++;
    }
  
  for ( i = 0; i < numberOfParameters; i++ )
    {
    if ( strcmp( name, tempNames[i] ) == 0 )
      {
      if ( this->TransformParameters[i] != value )
        {
        value = this->TransformParameters[i];
        }
      }
    }

  return value;
}

//----------------------------------------------------------------------------
vtkAbstractTransform*
vtkImageRegistration::GetLastTransform()
{
  switch ( this->TransformType )
    {
    case VTK_TRANSFORM_CENTERED:
      {
      vtkCenteredTransform *transform = 
        dynamic_cast<vtkCenteredTransform*>(this->Transform);
      vtkAmoebaMinimizer *optimizer = 
        dynamic_cast<vtkAmoebaMinimizer*>(this->Optimizer);

      for(int i = 0; i < 10; i++)
        {
        this->TransformParameters[i] = optimizer->GetParameterValue(
          this->GetTransformParameterName(this->TransformType, i));
        }

      transform->SetCenter( this->TransformParameters[0],
                            this->TransformParameters[1],
                            this->TransformParameters[2] );
      transform->SetTranslation( this->TransformParameters[3],
                                 this->TransformParameters[4],
                                 this->TransformParameters[5] );
      transform->SetRotationAnglesYXZ( this->TransformParameters[6],
                                       this->TransformParameters[7],
                                       this->TransformParameters[8] );
      transform->SetIsotropicScale( this->TransformParameters[9] );
      transform->Update();
      }
    default:
      vtkErrorMacro(<< "GetLastTransform: Unknown TransformType");
    }     

  return static_cast<vtkAbstractTransform*>(this->Transform);
}

//--------------------------------------------------------------------------
// Template function 
void 
vtkImageRegistration::InitializeTransformAndParameters()
{
  if (this->Transform != NULL)
    {
    this->Transform->Delete();
    this->Transform = NULL;
    }

  switch (this->TransformType)
    {
    case VTK_TRANSFORM_CENTERED:
      {
      vtkCenteredTransform *transform;

      transform = vtkCenteredTransform::New();
      transform->SetCenter( this->TransformParameters[0],
                            this->TransformParameters[1],
                            this->TransformParameters[2] );
      transform->SetTranslation( this->TransformParameters[3],
                                 this->TransformParameters[4],
                                 this->TransformParameters[5] );
      transform->SetRotationAnglesYXZ( this->TransformParameters[6],
                                       this->TransformParameters[7],
                                       this->TransformParameters[8] );
      transform->SetIsotropicScale( this->TransformParameters[9] );
      this->Transform = transform;
      }
      break;
    default:
      vtkErrorMacro( << "incorrect transform type!" );
      break;
    }
}

//--------------------------------------------------------------------------
// Template function 
void 
vtkImageRegistration::InitializeInterpolator()
{
  if (this->Interpolator == NULL)
    {
    this->Interpolator = vtkImageReslice::New();
    }
  this->Interpolator->SetInput(this->GetMovingImage());
  this->Interpolator->SetInformationInput(this->GetFixedImage());
  this->Interpolator->SetResliceTransform(this->Transform);

  switch(this->InterpolatorType)
    {
    case VTK_INTERPOLATOR_NEAREST_NEIGHBOR:
      {
      this->Interpolator->SetInterpolationModeToNearestNeighbor();
      }
      break;
    
    case VTK_INTERPOLATOR_LINEAR:
      {
      this->Interpolator->SetInterpolationModeToLinear();
      }
      break;

    case VTK_INTERPOLATOR_CUBIC:
      {
      this->Interpolator->SetInterpolationModeToCubic();
      }
      break;
    default:
      vtkErrorMacro( << "incorrect interpolater type!" );
      break;
    }
  this->Interpolator->OptimizationOff();
}

//--------------------------------------------------------------------------
// 
void 
vtkImageRegistration::InitializeMetricAndParameters()
{
  if (this->Metric != NULL)
    {
    this->Metric->Delete();
    this->Metric = NULL;
    }

  switch (this->MetricType)
    {
    case VTK_METRIC_NORMALIZED_CROSS_CORRELATION:
      {
      vtkCalcCrossCorrelation* metric;

      metric = vtkCalcCrossCorrelation::New();
      metric->SetInput1(this->GetFixedImage());
      metric->SetInput2(this->Interpolator->GetOutput());
      this->Metric = metric;
      }
      break;

    case VTK_METRIC_NORMALIZED_MUTUAL_INFORMATION:
      {
      vtkImageMutualInformation* metric;

      metric = vtkImageMutualInformation::New();
      metric->SetInput1(this->GetFixedImage());
      metric->SetInput2(this->Interpolator->GetOutput());
      this->Metric = metric;
      }
      break;
    default:
      vtkErrorMacro( << "incorrect metric type!" );
      break;
    }

}

//--------------------------------------------------------------------------
// 
static void vtkEvaluateFunction(void * arg)
{
  double val = 0.0;
  vtkImageRegistration::RegistrationInfo* registrationInfo = 
    static_cast<vtkImageRegistration::RegistrationInfo*>(arg);
  vtkImageRegistration* registration = registrationInfo->ImageRegistration;
  vtkObject* Optimizer = registrationInfo->Optimizer;
  vtkAbstractTransform* Transform = registrationInfo->Transform;
  vtkObject* Metric = registrationInfo->Metric;
  vtkImageReslice* Interpolator = registrationInfo->Interpolator;
#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION < 4)
  std::vector<double>* transformParametersPointer 
    = registrationInfo->TransformParametersPointer;
#else
  vtkstd::vector<double>* transformParametersPointer 
    = registrationInfo->TransformParametersPointer;
#endif

  switch (registrationInfo->TransformType)
    {
    case VTK_TRANSFORM_CENTERED:
      {
      vtkCenteredTransform *transform = 
        dynamic_cast<vtkCenteredTransform*>(Transform);
      vtkAmoebaMinimizer *optimizer = 
        dynamic_cast<vtkAmoebaMinimizer*>(Optimizer);

      for(int i = 3; i < 10; i++)
        {
        (*transformParametersPointer)[i] = optimizer->GetParameterValue(registration->GetTransformParameterName(registrationInfo->TransformType, i));
        }
      /*
      cerr << "Translation from amoeba: " 
           << (*transformParametersPointer)[3] << " " 
           << (*transformParametersPointer)[4] << " " 
           << (*transformParametersPointer)[5] << " " 
           << (*transformParametersPointer)[6] << " " 
           << (*transformParametersPointer)[7] << " " 
           << (*transformParametersPointer)[8] << " " 
           << (*transformParametersPointer)[9] << " "
           << endl;
      */     
      transform->SetTranslation( (*transformParametersPointer)[3],
                                 (*transformParametersPointer)[4],
                                 (*transformParametersPointer)[5] );
      transform->SetRotationAnglesYXZ( (*transformParametersPointer)[6],
                                       (*transformParametersPointer)[7],
                                       (*transformParametersPointer)[8] );
      transform->SetIsotropicScale( (*transformParametersPointer)[9] );
      transform->Update();
      Interpolator->UpdateWholeExtent();
      
      if (registrationInfo->MetricType == 
          VTK_METRIC_NORMALIZED_MUTUAL_INFORMATION)
        {
        vtkImageMutualInformation *metric = 
          dynamic_cast<vtkImageMutualInformation*>(Metric);

        metric->UpdateWholeExtent();

        val = metric->GetNormalizedMI();
        }
      else if (registrationInfo->MetricType == 
               VTK_METRIC_NORMALIZED_CROSS_CORRELATION)
        {
        vtkCalcCrossCorrelation *metric = 
          dynamic_cast<vtkCalcCrossCorrelation*>(Metric);
        
        metric->Update();

        val = metric->GetCrossCorrelation();
        }

      optimizer->SetFunctionValue(-val);
      }
      break;
    }
}

//------------------------------------------------------------------
// 
void 
vtkImageRegistration::InitializeOptimizerAndParameters()
{
  if (this->Optimizer != NULL)
    {
    this->Optimizer->Delete();
    this->Optimizer = NULL;
    }

  switch (this->OptimizerType)
    {
    case VTK_OPTIMIZER_AMOEBA:
      {
      vtkAmoebaMinimizer* optimizer;

      optimizer = vtkAmoebaMinimizer::New();
 
      optimizer->Initialize();
      optimizer->SetMaxIterations(100); 
      
      RegistrationArgs.ImageRegistration = this;
      RegistrationArgs.Optimizer = optimizer;
      RegistrationArgs.Transform = this->Transform;
      RegistrationArgs.Metric = this->Metric;
      RegistrationArgs.Interpolator = this->Interpolator;
      RegistrationArgs.OptimizerType = this->OptimizerType;
      RegistrationArgs.MetricType = this->MetricType;
      RegistrationArgs.InterpolatorType = this->InterpolatorType;
      RegistrationArgs.TransformType = this->TransformType;
      RegistrationArgs.TransformParametersPointer = &this->TransformParameters;

      optimizer->SetFunction(&vtkEvaluateFunction, 
                             (void*)(&this->RegistrationArgs));

      switch (this->TransformType)
        {
        case VTK_TRANSFORM_CENTERED:
          {
          double ts = this->TransformParameters[10];
          double rs = this->TransformParameters[11];
          double ss = this->TransformParameters[12];
          for (int i=3; i<10; i++)
            {
            optimizer->SetParameterValue(
              this->GetTransformParameterName(this->TransformType, i), 
              this->TransformParameters[i]);
            }
          
          optimizer->SetParameterScale(
            this->GetTransformParameterName(this->TransformType, 3), 
            ts); // translation x
          optimizer->SetParameterScale(
            this->GetTransformParameterName(this->TransformType, 4), 
            ts); // translation y
          optimizer->SetParameterScale(
            this->GetTransformParameterName(this->TransformType, 5), 
            ts); // translation z
          optimizer->SetParameterScale(
            this->GetTransformParameterName(this->TransformType, 6), 
            rs); // rotation x
          optimizer->SetParameterScale(
            this->GetTransformParameterName(this->TransformType, 7), 
            rs); // rotation y
          optimizer->SetParameterScale(
            this->GetTransformParameterName(this->TransformType, 8), 
            rs); // rotation z
          optimizer->SetParameterScale(
            this->GetTransformParameterName(this->TransformType, 9), 
            ss); // scaling
          }
          break;
        default:
          break;
        }
      this->Optimizer = optimizer;
      }
      break;
    default:
      vtkErrorMacro( << "incorrect optimizer type!" );
      break;
    }
}

//--------------------------------------------------------------------------
// 
int 
vtkImageRegistration::ExecuteRegistration()
{
  switch (this->OptimizerType)
    {
    case VTK_OPTIMIZER_AMOEBA:
      {
      dynamic_cast<vtkAmoebaMinimizer*>(this->Optimizer)->Minimize();
      /*
      for (int i=0; i<10; i++)
        {
        dynamic_cast<vtkAmoebaMinimizer*>(this->Optimizer)->Iterate();
        this->UpdateProgress(i/50.0);
        }
      */
      }
      break;
    default:
      vtkErrorMacro( << "incorrect optimizer type!" );
      break;
    }

  return 0;
}

//--------------------------------------------------------------------------
// The main program to set up the registration framework and connect
// pieces together and do registration
int
vtkImageRegistration::UpdateRegistration()
{
  vtkImageData *fixedImage = this->GetFixedImage();
  vtkImageData *movingImage = this->GetMovingImage();

  if (fixedImage->GetMTime() > this->ExecuteTime || 
      movingImage->GetMTime() > this->ExecuteTime || 
      this->GetMTime() > this->ExecuteTime )
    {
    int ret =  this->Initialize();
    if( ret != 0 )
      {
      return 1;
      }

    this->InvokeEvent(vtkCommand::StartEvent,NULL);
  
    // reset Abort flag
    this->AbortExecute = 0;
    this->Progress = 0.0;
    this->ExecuteRegistration();
    this->ExecuteTime.Modified();

    if ( !this->AbortExecute )
      {
      this->UpdateProgress(1.0);
      }
    
    this->InvokeEvent(vtkCommand::EndEvent,NULL);
    }

  return 0;
}

//--------------------------------------------------------------------------
// this is not an initialization routine rather a rountine to check
// for correct data type and components
int  
vtkImageRegistration::Initialize()
{
  vtkImageData *fixedImage = this->GetFixedImage();
  vtkImageData *movingImage = this->GetMovingImage();

  if ( !fixedImage )
    {
    vtkErrorMacro(<<"FixedImage is not present");
    return 1;
    }

  fixedImage->Update();

  if ( !movingImage )
    {
    vtkErrorMacro(<<"MovingImage is not present");
    return 1;
    }

  movingImage->Update();

  if ( fixedImage->GetScalarType() != movingImage->GetScalarType() )
    {
    vtkErrorMacro(<<"FixedImage type does not match MovingImage type");
    return 1;
    }

  /*
  vtkImageStencilData *stencil = this->GetFixedImageStencil();

  if (stencil)
    {
    stencil->SetSpacing(fixedImage->GetSpacing());
    stencil->SetOrigin(fixedImage->GetOrigin());
    stencil->SetUpdateExtent(stencil->GetWholeExtent());
    stencil->Update();
    }
  */

  if ( this->InterpolatorType < 0 || 
       this->InterpolatorType >= VTK_NUMBER_OF_INTERPOLATORS)
    {
    vtkErrorMacro(<<"Interpolator is not present" );
    return 1;
    }

  if ( this->MetricType < 0 || this->MetricType >= VTK_NUMBER_OF_METRICS )
    {
    vtkErrorMacro(<<"Metric is not present" );
    return 1;
    }

  if ( this->OptimizerType < 0 || 
       this->OptimizerType >= VTK_NUMBER_OF_OPTIMIZERS )
    {
    vtkErrorMacro(<<"Optimizer is not present" );
    return 1;
    }

  if( this->TransformType < 0 || 
      this->TransformType >= VTK_NUMBER_OF_TRANSFORMS )
    {
    vtkErrorMacro(<<"Transform is not present");
    return 1;
    }

  this->InitializeTransformAndParameters();

  this->InitializeInterpolator();
  
  this->InitializeMetricAndParameters();
  
  this->InitializeOptimizerAndParameters();
  
  return 0;
}
