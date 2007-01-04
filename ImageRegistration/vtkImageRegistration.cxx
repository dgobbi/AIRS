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
#include "vtkMatrix4x4.h"
#include "vtkHomogeneousTransform.h"
#include "vtkMatrixToHomogeneousTransform.h"
#include "vtkImageReslice.h"
#include "vtkImageGaussianSmooth.h"
#include "vtkImageShiftScale.h"
#include "vtkImageAccumulate.h"
#include "vtkCommand.h"
#if (VTK_MAJOR_VERSION >= 5) 
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#endif

// 
#include "vtkImageProjection.h"
#include "vtkImageRangeCalculator.h"

// Optimizer header files
#include "vtkAmoebaMinimizer.h"

// Image metric header files
#include "vtkCalcCrossCorrelation.h"
#include "vtkImageMutualInformation.h"

// Transform header files
#include "vtkCenteredTransform.h"

#ifndef vtkFloatingPointType
#define vtkFloatingPointType vtkFloatingPointType
typedef float vtkFloatingPointType;
#endif

const char * vtkOptimizerNames[] = {
  "vtkAmoebaMinimizer", 
  0
};

const char * vtkAmoebaMinimizerParameterNames[] = {
  "MaxNumberOfIterations",
  "Tolerance",
  0
};

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

const char * vtkInterpolatorNames[] = {
  "vtkNearestNeighborInterpolator",
  "vtkLinearInterpolator",
  "vtkCubicInterpolator", 
  0
};

const char * vtkPreprocessorNames[] = {
  "vtkEpilepsyMIPreprocessor",
  "vtkEpilepsyNCPreprocessor", 
  0
};

const char * vtkEpilepsyMIPreprocessorParameterNames[] = {
  "Blurring",
  "GaussScaleX",
  "GaussScaleY",
  "GaussScaleZ",
  "LowBin",
  "HighBin",
  "DownSampleFactorX",
  "DownSampleFactorY", 
  "DownSampleFactorZ",
  "SourceBackgroundLevel",
  "TargetBackgroundLevel",
  0
};

const char * vtkEpilepsyNCPreprocessorParameterNames[] = {
  "Blurring",
  "GaussScaleX",
  "GaussScaleY",
  "GaussScaleZ",
  "LowBin",
  "HighBin",
  "DownSampleFactorX",
  "DownSampleFactorY", 
  "DownSampleFactorZ",
  "SourceBackgroundLevel",
  "TargetBackgroundLevel",
  0
};

//----------------------------------------------------------------------------
vtkImageRegistration* vtkImageRegistration::New()
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
  this->PreprocessorType          = VTK_TYPE_INVALID;

  this->Transform                 = NULL; 
  this->Metric                    = NULL; 
  this->Optimizer                 = NULL; 

  this->SourceBlur                = NULL;
  this->TargetBlur                = NULL;
  this->SourceAccumulate          = NULL;
  this->TargetAccumulate          = NULL;
  this->SourceRange               = NULL;
  this->TargetRange               = NULL;
  this->SourceReslice             = NULL;
  this->TargetReslice             = NULL;
  this->SourceRescale             = NULL;
  this->TargetRescale             = NULL;

  this->LastTransform             = vtkTransform::New();

#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION < 4)
  this->OptimizerParameters       = std::vector<double>();
  this->MetricParameters          = std::vector<double>();
  this->TransformParameters       = std::vector<double>();
  this->PreprocessorParameters    = std::vector<double>();
#else
  this->OptimizerParameters       = vtkstd::vector<double>();
  this->MetricParameters          = vtkstd::vector<double>();
  this->TransformParameters       = vtkstd::vector<double>();
  this->PreprocessorParameters    = vtkstd::vector<double>();
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

  this->SourceImageRescaleIntercept   = 0.0;
  this->SourceImageRescaleSlope       = 1.0;
  this->TargetImageRescaleIntercept   = 0.0;
  this->TargetImageRescaleSlope       = 1.0;

#if (VTK_MAJOR_VERSION >= 5) 
  // we have the image inputs and the optional stencil input
  this->SetNumberOfInputPorts(2);
#endif
}

//----------------------------------------------------------------------------
vtkImageRegistration::~vtkImageRegistration() 
{
  // delete vtk objects
  if ( this->LastTransform != NULL )
    {
    this->LastTransform->Delete();
    }
  if ( this->SourceBlur != NULL )
    {
    this->SourceBlur->Delete();
    }
  if ( this->TargetBlur != NULL )
    {
    this->TargetBlur->Delete();
    }
  if ( this->SourceAccumulate != NULL )
    {
    this->SourceAccumulate->Delete();
    }
  if ( this->TargetAccumulate != NULL )
    {
    this->TargetAccumulate->Delete();
    }
  if ( this->SourceRange != NULL )
    {
    this->SourceRange->Delete();
    }
  if ( this->TargetRange != NULL )
    {
    this->TargetRange->Delete();
    }
  if ( this->SourceRescale != NULL )
    {
    this->SourceRescale->Delete();
    }
  if ( this->TargetRescale != NULL )
    {
    this->TargetRescale->Delete();
    }
  if ( this->SourceReslice != NULL )
    {
    this->SourceReslice->Delete();
    }
  if ( this->TargetReslice != NULL )
    {
    this->TargetReslice->Delete();
    }
  if ( this->Optimizer != NULL )
    {
    this->Optimizer->Delete();
    }
  if ( this->Transform != NULL )
    {
    this->Transform->Delete();
    }
  if ( this->Metric != NULL )
    {
    this->Metric->Delete();
    }

}
  
//----------------------------------------------------------------------------
void vtkImageRegistration::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkProcessObject::PrintSelf(os,indent);
  
  os << indent << "Optimizer: " << this->OptimizerType << "\n";
  os << indent << "Metric: " << this->MetricType << "\n";
  os << indent << "Interpolator: " << this->InterpolatorType << "\n";
  os << indent << "Transform: " << this->TransformType << "\n";
  os << indent << "Preprocessor: " << this->PreprocessorType << "\n";
}  

//----------------------------------------------------------------------------
void vtkImageRegistration::SetFixedImage(vtkImageData *input)
{
#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION <= 4)
  this->vtkProcessObject::SetNthInput(0, input);
#else
  // Ask the superclass to connect the input.
  this->SetNthInputConnection(0, 0, (input ? input->GetProducerPort() : 0));
#endif
}

//----------------------------------------------------------------------------
vtkImageData* vtkImageRegistration::GetFixedImage()
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
void vtkImageRegistration::SetMovingImage(vtkImageData *input)
{
#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION <= 4)
  this->vtkProcessObject::SetNthInput(1, input);
#else
  // Ask the superclass to connect the input.
  this->SetNthInputConnection(0, 1, (input ? input->GetProducerPort() : 0));
#endif
}

//----------------------------------------------------------------------------
vtkImageData* vtkImageRegistration::GetMovingImage()
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
void vtkImageRegistration::SetFixedImageStencil(vtkImageStencilData *stencil)
{
#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION <= 4)
  this->vtkProcessObject::SetNthInput(2, stencil);
#else
  // if stencil is null, then set the input port to null
  this->SetNthInputConnection(1, 0, 
    (stencil ? stencil->GetProducerPort() : 0));
#endif
}

//----------------------------------------------------------------------------
vtkImageStencilData* vtkImageRegistration::GetFixedImageStencil()
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
  if (this->GetNumberOfInputConnections(1) < 1)
    {
    return NULL;
    }
  return vtkImageStencilData::SafeDownCast(
    this->GetExecutive()->GetInputData(1, 0));
#endif
}

//----------------------------------------------------------------------------
const char* vtkImageRegistration::GetOptimizerName(int i)
{
  if ( i < 0 || i > VTK_NUMBER_OF_OPTIMIZERS )
    {
    return NULL;
    }
 
  return vtkOptimizerNames[i];
}

//----------------------------------------------------------------------------
const char* vtkImageRegistration::GetOptimizerParameterName(int optimizer, 
                                                            int parameter)
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
void vtkImageRegistration::SetOptimizerParameter(const char * name, 
                                                 double value)
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
const char* vtkImageRegistration::GetMetricName(int i)
{
  if ( i < 0 || i > VTK_NUMBER_OF_METRICS )
    {
    return NULL;
    }

  return vtkMetricNames[i];
}

//----------------------------------------------------------------------------
const char* vtkImageRegistration::GetMetricParameterName(int metric, 
                                                         int parameter)
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
void vtkImageRegistration::SetMetricParameter(const char * name , double value)
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
const char* vtkImageRegistration::GetPreprocessorName(int i)
{
  if ( i < 0 || i > VTK_NUMBER_OF_PREPROCESSORS )
    {
    return NULL;
    }
 
  return vtkPreprocessorNames[i];
}

//----------------------------------------------------------------------------
const char* vtkImageRegistration::GetPreprocessorParameterName(int preprocessor, 
                                                               int parameter)
{
  int i, numberOfParameters;
  const char ** tempNames = NULL;

  if ( preprocessor < 0 || preprocessor > VTK_NUMBER_OF_PREPROCESSORS )
    {
    return NULL;
    }
  
  switch ( preprocessor )
    {
    case VTK_OPTIMIZER_AMOEBA:
      tempNames = vtkEpilepsyMIPreprocessorParameterNames;
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
void vtkImageRegistration::SetPreprocessorParameter(const char * name, 
                                                    double value)
{
  unsigned int i, numberOfParameters;
  const char ** tempNames = NULL;

  if ( this->PreprocessorType < 0 ||  this->PreprocessorType 
       >= VTK_NUMBER_OF_PREPROCESSORS )
    {
    vtkErrorMacro( << "SetPreprocessorParameter: set preprocessor type before set optimizer parameters");
    return ;
    }
  
  switch ( this->PreprocessorType )
    {
    case VTK_PREPROCESSOR_MI:
      tempNames = vtkEpilepsyMIPreprocessorParameterNames;
      break;
    case VTK_PREPROCESSOR_NC:
      tempNames = vtkEpilepsyNCPreprocessorParameterNames;
      break;
    }     

  numberOfParameters = 0;
  for ( i = 0; tempNames[i] != 0; i++ )
    {
    numberOfParameters++;
    }
  
  if ( this->PreprocessorParameters.size() != numberOfParameters ) 
    {
    this->PreprocessorParameters.resize( numberOfParameters );
    }
  
  for ( i = 0; i < numberOfParameters; i++ )
    {
    if ( strcmp( name, tempNames[i] ) == 0 )
      {
      if ( this->PreprocessorParameters[i] != value )
        {
        this->PreprocessorParameters[i] = value;
        this->Modified();
        }
      }
    } 
}

//----------------------------------------------------------------------------
const char* vtkImageRegistration::GetInterpolatorName(int i)
{
  if ( i < 0 || i > VTK_NUMBER_OF_INTERPOLATORS )
    {
    return NULL;
    }

  return vtkInterpolatorNames[i];
}

//----------------------------------------------------------------------------
const char* vtkImageRegistration::GetTransformName(int i)
{
  if ( i < 0 || i > VTK_NUMBER_OF_TRANSFORMS )
    {
    return NULL;
    }

  return vtkTransformNames[i];
}

const char * vtkImageRegistration::GetTransformParameterName(int transform, 
                                                             int parameter)
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
void vtkImageRegistration::SetTransformParameter(const char * name, 
                                                 double value)
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
double vtkImageRegistration::GetTransformParameter(const char * name)
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
vtkTransform* vtkImageRegistration::GetLastTransform()
{
  this->LastTransform->Identity();
  switch ( this->TransformType )
    {
    case VTK_TRANSFORM_CENTERED:
      {
      vtkCenteredTransform *transform = 
        dynamic_cast<vtkCenteredTransform*>(this->Transform);
      vtkAmoebaMinimizer *optimizer = 
        dynamic_cast<vtkAmoebaMinimizer*>(this->Optimizer);

      for(int i = 3; i < 10; i++)
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
      this->LastTransform->SetMatrix(transform->GetMatrix());
      }
      break;
    default:
      vtkErrorMacro(<< "GetLastTransform: Unknown TransformType");
    }     

  return this->LastTransform;
}

//--------------------------------------------------------------------------
void vtkImageRegistration::InitializeTransform()
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
      transform->Update();
      this->Transform = transform;
      }
      break;
    default:
      vtkErrorMacro( << "incorrect transform type!" );
      break;
    }
}

//--------------------------------------------------------------------------
void vtkImageRegistration::InitializePreprocessor(void)
{
  int sourceDim[3], targetDim[3], downsampleFactor[3], resolution[3];
  vtkFloatingPointType sourceSpacing[3], targetSpacing[3];
  vtkFloatingPointType targetOrigin[3];
  vtkFloatingPointType targetScale[3];

  if (this->SourceBlur != NULL)
    {
    this->SourceBlur->Delete();
    }
  this->SourceBlur = vtkImageGaussianSmooth::New();
  
  if (this->TargetBlur != NULL)
    {
    this->TargetBlur->Delete();
    }
  this->TargetBlur = vtkImageGaussianSmooth::New();

  if (this->SourceAccumulate != NULL)
    {
    this->SourceAccumulate->Delete();
    }
  this->SourceAccumulate = vtkImageAccumulate::New();
  
  if (this->TargetAccumulate != NULL)
    {
    this->TargetAccumulate->Delete();
    }
  this->TargetAccumulate = vtkImageAccumulate::New();

  if (this->SourceRange != NULL)
    {
    this->SourceRange->Delete();
    }
  this->SourceRange = vtkImageRangeCalculator::New();
  
  if (this->TargetRange != NULL)
    {
    this->TargetRange->Delete();
    }
  this->TargetRange = vtkImageRangeCalculator::New();

  if (this->SourceRescale != NULL)
    {
    this->SourceRescale->Delete();
    }
  this->SourceRescale = vtkImageShiftScale::New();

  if (this->TargetRescale != NULL)
    {
    this->TargetRescale->Delete();
    }
  this->TargetRescale = vtkImageShiftScale::New();
  
  if (this->SourceReslice != NULL)
    {
    this->SourceReslice->Delete();
    }
  this->SourceReslice = vtkImageReslice::New();

  if (this->TargetReslice != NULL)
    {
    this->TargetReslice->Delete();
    }
  this->TargetReslice = vtkImageReslice::New();

  // Connect pipeline
  if (this->PreprocessorParameters[0] == 0)
    {
    this->SourceReslice->SetInput(this->GetMovingImage());
    this->TargetReslice->SetInput(this->GetFixedImage());
    }
  else
    {
    this->SourceBlur->SetInput(this->GetMovingImage());
    this->TargetBlur->SetInput(this->GetFixedImage());

    this->SourceReslice->SetInput(this->SourceBlur->GetOutput());
    this->TargetReslice->SetInput(this->TargetBlur->GetOutput());
    }

  this->SourceRescale->SetInput(this->SourceReslice->GetOutput());
  this->TargetRescale->SetInput(this->TargetReslice->GetOutput());

  // Set up reslice parameters
  switch(this->InterpolatorType)
    {
    case VTK_INTERPOLATOR_NEAREST_NEIGHBOR:
      {
      this->SourceReslice->SetInterpolationModeToNearestNeighbor();
      this->TargetReslice->SetInterpolationModeToNearestNeighbor();
      }
      break;
    
    case VTK_INTERPOLATOR_LINEAR:
      {
      this->SourceReslice->SetInterpolationModeToLinear();
      this->TargetReslice->SetInterpolationModeToLinear();
      }
      break;

    case VTK_INTERPOLATOR_CUBIC:
      {
      this->SourceReslice->SetInterpolationModeToCubic();
      this->TargetReslice->SetInterpolationModeToCubic();
      }
      break;
    default:
      vtkErrorMacro( << "incorrect interpolater type!" );
      break;
    }
  this->SourceReslice->OptimizationOff();
  this->TargetReslice->OptimizationOff();
  this->SourceReslice->SetResliceTransform(this->Transform);

  this->GetMovingImage()->Update();
  this->GetMovingImage()->GetDimensions(sourceDim);
  this->GetFixedImage()->Update();
  this->GetFixedImage()->GetDimensions(targetDim);

  this->GetMovingImage()->GetSpacing(sourceSpacing);
  this->GetFixedImage()->GetSpacing(targetSpacing);
  
  this->GetFixedImage()->GetOrigin(targetOrigin);

  downsampleFactor[0] = (int)(this->PreprocessorParameters[6]);
  downsampleFactor[1] = (int)(this->PreprocessorParameters[7]);
  downsampleFactor[2] = (int)(this->PreprocessorParameters[8]);
  resolution[0] = targetDim[0]/downsampleFactor[0];
  resolution[1] = targetDim[1]/downsampleFactor[1];
  resolution[2] = targetDim[2]/downsampleFactor[2];

  targetScale[0] = (targetDim[0]-1.0)/(resolution[0]-1.0);
  targetScale[1] = (targetDim[1]-1.0)/(resolution[1]-1.0);
  targetScale[2] = (targetDim[2]-1.0)/(resolution[2]-1.0);

  this->SourceReslice->SetOutputExtent(0, resolution[0]-1,
                                       0, resolution[1]-1,
                                       0, resolution[2]-1);
  this->TargetReslice->SetOutputExtent(0, resolution[0]-1,
                                       0, resolution[1]-1,
                                       0, resolution[2]-1);
  this->SourceReslice->SetOutputSpacing(targetSpacing[0]*targetScale[0],
                                        targetSpacing[1]*targetScale[1],
                                        targetSpacing[2]*targetScale[2]);
  this->TargetReslice->SetOutputSpacing(targetSpacing[0]*targetScale[0],
                                        targetSpacing[1]*targetScale[1],
                                        targetSpacing[2]*targetScale[2]);

  this->SourceReslice->SetOutputOrigin(targetOrigin);
  this->TargetReslice->SetOutputOrigin(targetOrigin);

//   vtkImageStencilData *stencil = this->GetFixedImageStencil();
//   if (stencil)
//     {
//     this->SourceReslice->SetStencil(stencil);
//     this->TargetReslice->SetStencil(stencil);
//     }

  // compute the range
  vtkFloatingPointType sourceRange[2], targetRange[2];
  double sourceFractionRange[2], targetFractionRange[2];
  double sourceUpperValue, sourceLowerValue;
  double targetUpperValue, targetLowerValue;

  this->GetMovingImage()->GetScalarRange(sourceRange);

  this->SourceAccumulate->SetInput(this->GetMovingImage());
  this->SourceAccumulate->SetComponentSpacing(1.0, 0, 0);
  this->SourceAccumulate->SetComponentOrigin(sourceRange[0], 0, 0);
  this->SourceAccumulate->SetComponentExtent(0, int(0.5+(sourceRange[1]-sourceRange[0])), 0, 0, 0, 0);

  this->SourceRange->SetInput(this->SourceAccumulate->GetOutput());
  this->SourceRange->SetAreaFractionRange(0.0, 0.99);
  this->SourceRange->Calculate();
  this->SourceRange->GetDataRange(sourceFractionRange);

  
  sourceLowerValue = floor(this->PreprocessorParameters[9]);
  sourceUpperValue = sourceFractionRange[1];

  this->GetFixedImage()->GetScalarRange(targetRange);

  this->TargetAccumulate->SetInput(this->GetFixedImage());
  this->TargetAccumulate->SetComponentSpacing(1.0, 0, 0);
  this->TargetAccumulate->SetComponentOrigin(targetRange[0], 0, 0);
  this->TargetAccumulate->SetComponentExtent(0, int(0.5+(targetRange[1]-targetRange[0])), 0, 0, 0, 0);

  this->TargetRange->SetInput(this->TargetAccumulate->GetOutput());
  this->TargetRange->SetAreaFractionRange(0.0, 0.99);
  this->TargetRange->Calculate();
  this->TargetRange->GetDataRange(targetFractionRange);

  targetLowerValue = floor(this->PreprocessorParameters[10]);
  targetUpperValue = targetFractionRange[1];

  // Set up rescale parameters
  this->SourceRescale->SetOutputScalarTypeToUnsignedChar();
  this->TargetRescale->SetOutputScalarTypeToUnsignedChar();
  this->SourceRescale->ClampOverflowOn();
  this->TargetRescale->ClampOverflowOn();

  int lowBin = (int)(this->PreprocessorParameters[4]);
  int highBin = (int)(this->PreprocessorParameters[5]);

  this->SourceRescale->SetShift((lowBin+0.5)*(sourceUpperValue-sourceLowerValue)/(highBin-lowBin) - sourceLowerValue);
  this->SourceRescale->SetScale((highBin-lowBin)/(sourceUpperValue-sourceLowerValue));
  this->TargetRescale->SetShift((lowBin+0.5)*(targetUpperValue-targetLowerValue)/(highBin-lowBin) - targetLowerValue);
  this->TargetRescale->SetScale((highBin-lowBin)/(targetUpperValue-targetLowerValue));

  switch(this->PreprocessorType)
    {
    case VTK_PREPROCESSOR_MI:
      {
      this->SourceBlur->SetStandardDeviations(this->PreprocessorParameters[1], 
                                              this->PreprocessorParameters[2],
                                              this->PreprocessorParameters[3]);
      this->TargetBlur->SetStandardDeviations(this->PreprocessorParameters[1], 
                                              this->PreprocessorParameters[2],
                                              this->PreprocessorParameters[3]);
      }
      break;
    case VTK_PREPROCESSOR_NC:
      {
      this->SourceBlur->SetStandardDeviations(this->PreprocessorParameters[1], 
                                              this->PreprocessorParameters[2],
                                              this->PreprocessorParameters[3]);
      this->TargetBlur->SetStandardDeviations(this->PreprocessorParameters[1], 
                                              this->PreprocessorParameters[2],
                                              this->PreprocessorParameters[3]);
      }
      break;
    }

}
  
//--------------------------------------------------------------------------
void vtkImageRegistration::InitializeMetric(void)
{
  if (this->Metric != NULL)
    {
    this->Metric->Delete();
    this->Metric = NULL;
    }
  
  vtkImageStencilData *stencil = this->GetFixedImageStencil();

  switch (this->MetricType)
    {
    case VTK_METRIC_NORMALIZED_CROSS_CORRELATION:
      {
      vtkCalcCrossCorrelation* metric;

      metric = vtkCalcCrossCorrelation::New();
      metric->SetInput1(this->SourceRescale->GetOutput());
      metric->SetInput2(this->TargetRescale->GetOutput());
      if (stencil)
        {
        metric->SetStencil(stencil);
        }
      this->Metric = metric;
      }
      break;

    case VTK_METRIC_NORMALIZED_MUTUAL_INFORMATION:
      {
      vtkImageMutualInformation* metric;

      metric = vtkImageMutualInformation::New();
      metric->SetInput1(this->SourceRescale->GetOutput());
      metric->SetInput2(this->TargetRescale->GetOutput());
      metric->SetImageAComponentExtent(0, (int)(this->MetricParameters[0]+1));
      metric->SetImageBComponentExtent(0, (int)(this->MetricParameters[1]+1));
      
      if (stencil)
        {
        metric->SetStencil(stencil);
        }
      this->Metric = metric;
      }
      break;
    default:
      vtkErrorMacro( << "incorrect metric type!" );
      break;
    }
}

//--------------------------------------------------------------------------
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
      
//        cerr << "Translation from amoeba: " 
//             << (*transformParametersPointer)[3] << " " 
//             << (*transformParametersPointer)[4] << " " 
//             << (*transformParametersPointer)[5] << " " 
//             << (*transformParametersPointer)[6] << " " 
//             << (*transformParametersPointer)[7] << " " 
//             << (*transformParametersPointer)[8] << " " 
//             << (*transformParametersPointer)[9] << " "
//             << endl;
           
      transform->SetTranslation( (*transformParametersPointer)[3],
                                 (*transformParametersPointer)[4],
                                 (*transformParametersPointer)[5] );
      transform->SetRotationAnglesYXZ( (*transformParametersPointer)[6],
                                       (*transformParametersPointer)[7],
                                       (*transformParametersPointer)[8] );
      transform->SetIsotropicScale( (*transformParametersPointer)[9] );
      transform->Update();
      //Interpolator->UpdateWholeExtent();
      
      if (registrationInfo->MetricType == 
          VTK_METRIC_NORMALIZED_MUTUAL_INFORMATION)
        {
        vtkImageMutualInformation *metric = 
          dynamic_cast<vtkImageMutualInformation*>(Metric);

        metric->Update();

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
//       cerr << "metric val: " << val << endl;
      optimizer->SetFunctionValue(-val);
      }
      break;
    }
}

//------------------------------------------------------------------
void vtkImageRegistration::InitializeOptimizer()
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
      double tolerance;

      optimizer = vtkAmoebaMinimizer::New();
 
      optimizer->Initialize();
      // optimizer->SetMaxIterations((int)(this->OptimizerParameters[0])); 
      tolerance = this->OptimizerParameters[1];
      if (tolerance < 0.01 && 
	  tolerance > 0.000001)
	{
	optimizer->SetTolerance(tolerance);
	}

      RegistrationArgs.ImageRegistration = this;
      RegistrationArgs.Optimizer = optimizer;
      RegistrationArgs.Transform = this->Transform;
      RegistrationArgs.Metric = this->Metric;
      RegistrationArgs.Interpolator = this->SourceReslice;
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
int vtkImageRegistration::ExecuteRegistration()
{
  switch (this->OptimizerType)
    {
    case VTK_OPTIMIZER_AMOEBA:
      {
//       dynamic_cast<vtkAmoebaMinimizer*>(this->Optimizer)->Minimize();
      
      int iterNo = (int)(this->OptimizerParameters[0]);
      double tolerance = this->OptimizerParameters[1];
      if (iterNo != 0)
        {
        for (int i=0; i<iterNo; i++)
          {
          dynamic_cast<vtkAmoebaMinimizer*>(this->Optimizer)->Iterate();
          this->UpdateProgress(1.0*i/iterNo);
          }
        }
      else 
        {
        dynamic_cast<vtkAmoebaMinimizer*>(this->Optimizer)->Minimize();
        }
      }
      break;
    default:
      vtkErrorMacro( << "incorrect optimizer type!" );
      break;
    }

  return 0;
}

//--------------------------------------------------------------------------
int vtkImageRegistration::UpdateRegistration()
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
    
    //
//     this->SourceReslice->UpdateWholeExtent();
//     this->TargetReslice->UpdateWholeExtent();
    
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
int vtkImageRegistration::Initialize()
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

  vtkImageStencilData *stencil = this->GetFixedImageStencil();

  if (stencil)
    {
    stencil->SetSpacing(fixedImage->GetSpacing());
    stencil->SetOrigin(fixedImage->GetOrigin());
    stencil->SetUpdateExtent(fixedImage->GetWholeExtent());
    stencil->Update();
    }
 
  if ( this->InterpolatorType < 0 || 
       this->InterpolatorType >= VTK_NUMBER_OF_INTERPOLATORS)
    {
    vtkErrorMacro(<<"Interpolator is not present!" );
    return 1;
    }

  if ( this->MetricType < 0 || this->MetricType >= VTK_NUMBER_OF_METRICS )
    {
    vtkErrorMacro(<<"Metric is not present!" );
    return 1;
    }

  if ( this->PreprocessorType < 0 || 
       this->PreprocessorType >= VTK_NUMBER_OF_PREPROCESSORS)
    {
    vtkErrorMacro(<<"Preprocessor is not present!" );
    return 1;
    }

  if ( this->OptimizerType < 0 || 
       this->OptimizerType >= VTK_NUMBER_OF_OPTIMIZERS )
    {
    vtkErrorMacro(<<"Optimizer is not present!" );
    return 1;
    }

  if( this->TransformType < 0 || 
      this->TransformType >= VTK_NUMBER_OF_TRANSFORMS )
    {
    vtkErrorMacro(<<"Transform is not present!");
    return 1;
    }

  this->InitializeTransform();

  this->InitializePreprocessor();

  this->InitializeMetric();

  this->InitializeOptimizer();
  
  return 0;
}

#if (VTK_MAJOR_VERSION >= 5) 
//----------------------------------------------------------------------------
int vtkImageRegistration::FillInputPortInformation(int port, 
                                                   vtkInformation* info)
{
  if (port == 0)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    }
  if (port == 1)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageStencilData");
    // the stencil input is optional
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }
  return 1;
}
#endif

