/*=========================================================================

  Program:   Atamai Classes for VTK
  Module:    vtkITKImageSegmentation3D.cxx
  Creator:   Piali Das <pdas@atamai.com>

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

#include "vtkITKImageSegmentation3D.h"

// VTK header files
#include "vtkAlgorithmOutput.h"
#include "vtkDataObject.h"
#include "vtkMath.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPointSet.h"
#include "vtkTransform.h"
#include "vtkTimerLog.h"
#include "vtkExecutive.h"
#include "vtkObjectFactory.h"
#include "vtkImageData.h"
#include "vtkImageImport.h"
#include "vtkImageExport.h"
#include "vtkImageToImageStencil.h"
#include "vtkImageStencilData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkCommand.h"

// Headers for Watershed
#include "itkGradientMagnitudeImageFilter.h"
#include "itkWatershedImageFilter.h"

// Headers for Canny Segmentation Level Set Segmentation
#include "itkCannySegmentationLevelSetImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"

#include "itkFastMarchingImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkSigmoidImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkZeroCrossingImageFilter.h"

// knowledge based
#include "math.h"
#include "itkVector.h"
#include "itkBayesianClassifierImageFilter.h"
#include "itkVectorImage.h"
#include "itkBayesianClassifierInitializationImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

// itk image import/export files
#include "itkImage.h"
#include "itkVTKImageImport.h"
#include "itkVTKImageExport.h"
#include "ConnectVTKITK.h"
#include "itkCastImageFilter.h"

#include "itkCommand.h"


//--------------------------------------------------------------------------
//  NAMES OF THE ALGORITHMS
//--------------------------------------------------------------------------
const char *vtkAlgorithmNames[] = {
  "vtkCannySegmentationLevelSet",
  "vtkFastMarchLevelSet",
  "vtkWaterShed",
  "vtkKnowledgeBased",
   0
};

//--------------------------------------------------------------------------
//  Names Of the  Parameters for respective Segmentation Algorithms

const char *vtkCannySegmentationLevelSetParameterName[] = {
  "Threshold",
  "Variance",
  "CurvatureScaling",
  "PropagationScaling",
  "AdvectionScaling",
  "MaximumRMSError",
  "NumberOfIterationsofLevelSet",
  "IsoSurfaceValue",
  "NumberOfIterationsforDiffusion", //  Gradient Anisotropic Diffusion Filter
  "TimeStep",
  "ConductanceParameter",
  "LowerThreshold",
  "UpperThreshold",
  "OutsideValue",
  "InsideValue",
   0
};

const char *vtkFastMarchingSegmentationLevelSetParameterName[] = {
  "TimeStep",        // Diffusion filter
  "NumberOfIterationsOfDiffusion",
  "ConductanceParameter",
  "Sigma",           // Gradient filter
  "OutputMinimum",   // Sigmoid filter
  "OutputMaximum",
  "Alpha",
  "Beta",
  "StoppingValue",   // FastMarching Parameters
  "UpperThreshold",  // Threshold filter
  "LowerThreshold",
  "OutsideValue",
  "InsideValue",
  0
};

const char *vtkWatershedSegmentationParameterName[] = {
  "Level",
  "Thredhold",
  "TimeStep",
  "ConductanceParameter",
  "NumberOfIterations",
  0
};

const char *vtkKnowledgeBasedSegmentationParameterName[] = {
  "Conductance",
  "TimeStep",
  "NumberOfIterations",
  "OutputMinimum",
  "OutputMaximum",
  "NumberOfClasses",
  0
};

//----------------------------------------------------------------------------
// Implements the New method for the class
//----------------------------------------------------------------------------
vtkITKImageSegmentation3D* vtkITKImageSegmentation3D::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret =
    vtkObjectFactory::CreateInstance("vtkITKImageSegmentation3D");

  if(ret)
    {
    return (vtkITKImageSegmentation3D*)ret;
    }

  // If the factory was unable to create the object, then create it here.
  return new vtkITKImageSegmentation3D;
}


//----------------------------------------------------------------------------
//CONSTRUCTOR
//----------------------------------------------------------------------------
vtkITKImageSegmentation3D::vtkITKImageSegmentation3D()
{
  this->Algorithm = VTK_SEGMENTATION_ALGORITHM_INVALID;

  this->Parameters = itk::Array< double >(1);
  this->DefaultParameterValues = itk::Array< double >(1);
  this->VTKImageImporter = vtkImageImport::New();
  this->ImageToStencil = NULL;
  this->Initialized  = 0;
  this->FeatureImageDimension = 0;
  this->LabelImageDimension =0;
  this->FeatureImageDataType = 0;
  this->StencilThreshold = 1.0;
  this->ParameterChanged = 0;
  this->SetNumberOfInputPorts(3);
  this->SetNumberOfOutputPorts(2);
}

//----------------------------------------------------------------------------
//DESTRUCTOR
//----------------------------------------------------------------------------
vtkITKImageSegmentation3D::~vtkITKImageSegmentation3D()
{
  // delete itk objects
  if (this->ITKImageImporter != NULL)
    {
    this->ITKImageImporter->UnRegister();
    }

  if (this->ITKLabelImageImporter != NULL)
    {
    this->ITKLabelImageImporter->UnRegister();
    }

  for (int i = 0; i < static_cast<int>(Filters.size()); i++)
    {
    if (this->Filters.at(i) != NULL)
      {
      this->Filters.at(i)->UnRegister();
      }
    }

  if (this->VTKImageImporter)
    {
    this->VTKImageImporter->Delete();
    }

  if (this->ImageToStencil)
    {
    this->ImageToStencil->Delete();
    }
}

//---------------------------------------------------------------------------
 void vtkITKImageSegmentation3D::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "Algorithm: " << this->Algorithm << "\n";
}


//----------------------------------------------------------------------------
 void vtkITKImageSegmentation3D::SetInput( vtkImageData *input)
{
  if (input)
    {
    this->SetInputConnection(0, input->GetProducerPort());
    }
  else
    {
    this->SetInputConnection(0,0);
    }
}
//----------------------------------------------------------------------------
vtkImageData* vtkITKImageSegmentation3D::GetInput()
{
 return vtkImageData::SafeDownCast( this->GetExecutive()->GetInputData(0, 0));
}

//----------------------------------------------------------------------------
void vtkITKImageSegmentation3D::SetFeatureImage( vtkImageData *input)
{
  this->SetInput(input);
}

//----------------------------------------------------------------------------
vtkImageData *vtkITKImageSegmentation3D::GetFeatureImage()
{
  return vtkImageData::SafeDownCast(this->GetExecutive()->GetInputData(0, 0));
}

//----------------------------------------------------------------------------
void vtkITKImageSegmentation3D::SetLabelImage(vtkImageData *input)
{
  if (input)
    {
    this->SetInputConnection(1, input->GetProducerPort());
    }
  else
    {
    this->SetInputConnection(1,0);
    }
}
//---------------------------------------------------------------------------
vtkImageData *vtkITKImageSegmentation3D::GetLabelImage()
{
  return vtkImageData::SafeDownCast(this->GetExecutive()->GetInputData(1, 0));
}

//----------------------------------------------------------------------------
void vtkITKImageSegmentation3D::SetSeedPoints(vtkPointSet *input)
{
  if (input)
    {
    this->SetInputConnection(2, input->GetProducerPort());
    }
  else
    {
    this->SetInputConnection(2,0);
    }
}

//---------------------------------------------------------------------------
// returns the seed points
//---------------------------------------------------------------------------
vtkPointSet *vtkITKImageSegmentation3D::GetSeedPoints()
{
  return vtkPointSet::SafeDownCast(this->GetExecutive()->GetInputData(2, 0));
}

//----------------------------------------------------------------------------
int  vtkITKImageSegmentation3D::FillInputPortInformation(
  int port, vtkInformation* info)
{
  if (port == 0)
    {
    // The feature image input
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    }
  else if (port == 1)
    {
    // The label image input
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }
  else if (port == 2)
    {
    // seed point input
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }
  return 1;
}

//----------------------------------------------------------------------------
int  vtkITKImageSegmentation3D::FillOutputPortInformation(
  int port, vtkInformation* info)
{
  if (port == 0)
    {
    // The feature image input
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    }
  else if (port == 1)
    {
    // The label image input
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageStencilData");
    }
  return 1;
}

//----------------------------------------------------------------------------
// Gets Output
//----------------------------------------------------------------------------
vtkImageData *vtkITKImageSegmentation3D::GetOutput()
{
  return vtkImageData::SafeDownCast(this->GetOutputDataObject(0));
}

//----------------------------------------------------------------------------
vtkImageStencilData *vtkITKImageSegmentation3D::GetStencilOutput()
{
  if (this->ImageToStencil == NULL)
    {
    this->ImageToStencil = vtkImageToImageStencil::New();
    }
  vtkImageToImageStencil *imageToStencil = this->ImageToStencil;

  imageToStencil->SetInput(this->GetOutput());
  imageToStencil->SetLowerThreshold(0);
  imageToStencil->SetUpperThreshold(500.0);
  imageToStencil->Update();
  vtkImageStencilData *stencil =
    vtkImageStencilData::SafeDownCast(imageToStencil->GetOutput());

  return stencil;
}

//---------------------------------------------------------------------------
int vtkITKImageSegmentation3D::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inInfo,
  vtkInformationVector *outInfoVec)
{
  vtkImageData *input1 = vtkImageData::SafeDownCast(
    inInfo[0]->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));
  this->FeatureImageDataType = input1->GetScalarType();

  if(this->Algorithm != VTK_SEGMENTATION_ALGORITHM_INVALID)
    {
    if(this->GetNeedsLabelImage(this->GetAlgorithm()))
      {
      vtkImageData *input2 = vtkImageData::SafeDownCast(
        inInfo[1]->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));
      this->LabelImageDimension = input2->GetDataDimension();
      this->FeatureImageDimension = input1->GetDataDimension();
      if(this->FeatureImageDimension != this->FeatureImageDimension)
        {
        vtkErrorMacro(<<"Label Image Dimension does not match"
                      <<"Featur Image Dimension.");
        return 0;
        }
      }
    }

  // set default parameter values based on the feature image &
  // segmentation algorithm
  if(!this->SetDefaultParameterValues())
    {
    vtkErrorMacro(<<"vtkITKImageSegmentation3D::RequestInformation -"
                  <<"Algorithm or the feature image is not set properly.");
    }
  vtkInformation* outInfo1 = outInfoVec->GetInformationObject(0);
  vtkInformation* outInfo2 = outInfoVec->GetInformationObject(1);
  double dataSpacing[3], dataOrigin[3];
  int dataExtent[6];
  input1->GetSpacing(dataSpacing);
  input1->GetOrigin( dataOrigin);
  input1->GetWholeExtent(dataExtent);

  outInfo1->Set(vtkDataObject::SPACING(),dataSpacing,3);
  outInfo1->Set(vtkDataObject::ORIGIN(),dataOrigin,3);

  outInfo1->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
               dataExtent,6);
  outInfo2->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
                dataExtent,6);
  vtkDataObject::SetPointDataActiveScalarInfo(outInfo1,
                                              SEGMENTATION_OUTPUT_TYPE, 1);
  return 1;
}
//-----------------------------------------------------------------------------
int vtkITKImageSegmentation3D::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inInfoVec),
  vtkInformationVector *vtkNotUsed(outInfoVec))
{
  this->SetVTKITKPipelineConnection();

  return 1;
}

//-----------------------------------------------------------------------------
// Sets the Algorithm
//-----------------------------------------------------------------------------
void vtkITKImageSegmentation3D::SetAlgorithm(int algorithm)
{
  this->Algorithm = algorithm;
}

//-----------------------------------------------------------------------------
int vtkITKImageSegmentation3D::GetAlgorithm()
{
 return this->Algorithm;
}

//----------------------------------------------------------------------------
// returns the name of the algorithm
//----------------------------------------------------------------------------
const char *vtkITKImageSegmentation3D::GetAlgorithmName(int i)
{
  if (i < 0 || i >= VTK_NUMBER_OF_ALGORITHMS)
    {
    return NULL;
    }

  return vtkAlgorithmNames[i];
}

//----------------------------------------------------------------------------
const char* vtkITKImageSegmentation3D::GetTechniqueName(int techniqueID)
{
  return vtkAlgorithmNames[techniqueID];
}

//----------------------------------------------------------------------------
// Returns 1 if Label Image is required for a particular segmentation
// algorithm otherwise returns 0
int vtkITKImageSegmentation3D::GetNeedsLabelImage(int algorithm)
{
  int needLabelImage = -1;

  switch (algorithm)
    {
    case VTK_SEGMENTATION_ALGORITHM_INVALID:
      needLabelImage = 1;
      break;
    case VTK_LEVEL_SET_CANNY_EDGE:
      needLabelImage = 1;
      break;
    case VTK_LEVEL_SET_FAST_MARCHING:
      needLabelImage = 0;
      break;
    case VTK_WATERSHED:
      needLabelImage = 0;
      break;
    case VTK_KNOWLEDGEBASED:
      needLabelImage = 0;
      break;
    }

  return needLabelImage;
}

//----------------------------------------------------------------------------
// Returns 1 if Seeds are required for particular segmentation
// algorithm, otherwise returns 0
//----------------------------------------------------------------------------
int vtkITKImageSegmentation3D::GetNeedsSeedPoints(int algorithm)
{
  int needSeed = -1;

  switch (algorithm)
    {
    case VTK_SEGMENTATION_ALGORITHM_INVALID:
      needSeed = 0;
      break;
    case VTK_LEVEL_SET_CANNY_EDGE:
      needSeed = 0;
      break;
    case VTK_LEVEL_SET_FAST_MARCHING:
      needSeed = 1;
      break;
    case VTK_WATERSHED:
      needSeed = 0;
      break;
    case  VTK_KNOWLEDGEBASED:
      needSeed = 0;
      break;
    }

  return needSeed;
}
//-----------------------------------------------------------------------------
// returns number of inputs required by the algorithm
//-----------------------------------------------------------------------------
int vtkITKImageSegmentation3D::GetNumberOfInputs()
{
  int numberOfInputs = -1;

  switch (this->Algorithm)
    {
    case VTK_SEGMENTATION_ALGORITHM_INVALID:
      numberOfInputs = 0;
      break;
    case VTK_LEVEL_SET_CANNY_EDGE:
      numberOfInputs = 3;
      break;
    case VTK_LEVEL_SET_FAST_MARCHING:
      numberOfInputs = 2;
      break;
    case VTK_WATERSHED:
      numberOfInputs = 1;
      break;
    case  VTK_KNOWLEDGEBASED:
      numberOfInputs = 1;
      break;
    }

  return numberOfInputs;
}
//---------------------------------------------------------------------------
int vtkITKImageSegmentation3D::SetDefaultParameterValues()
{
  double min = 0.0;
  double max  = 0.0;
  int dim = 3;
  double defaultTimeStep = 0.0625;
  double defaultIsoSurface = 0.0;
  size_t numberOfParameters = 1;

  if (this->Algorithm == VTK_SEGMENTATION_ALGORITHM_INVALID)
    {
    return 0;
    }

    numberOfParameters = this->GetNumberOfParameters(this->Algorithm);
    vtkImageData *featureImage = vtkImageData::SafeDownCast(this->GetInput());
    vtkImageData *labelImage = vtkImageData::SafeDownCast(this->GetLabelImage());

    if (labelImage != NULL)
      {
      double  r[2];
      labelImage->GetScalarRange(r);
      defaultIsoSurface = (r[1]-r[0])/2;
      }

    if (featureImage != NULL)
      {
      min = 0.0;
      max = 255.0;
      dim = featureImage->GetDataDimension();
      }

    if (dim == 2)
      {
      defaultTimeStep = 0.125;
      }

    switch (this->Algorithm)
      {
      case VTK_LEVEL_SET_CANNY_EDGE:
        {
        if (featureImage != NULL)
          {
          if (this->DefaultParameterValues.Size() != numberOfParameters)
            {
            this->DefaultParameterValues.SetSize(numberOfParameters);
            }
          this->DefaultParameterValues[0] = 6.0;  //"Threshold",
          this->DefaultParameterValues[1] = 1.0;  //"Variance",
          this->DefaultParameterValues[2] = 1.0;  //"CurvatureScaling",
          this->DefaultParameterValues[3] = 0.0;  //"PropagationScaling",
          this->DefaultParameterValues[4] = 0.0;  //"AdvectionScaling",
          this->DefaultParameterValues[5] = 0.01;  //"MaximumRMSError",
          this->DefaultParameterValues[6] = 0.0;  //"NumberOfIterationsofLevelSet",
          this->DefaultParameterValues[7] = defaultIsoSurface;//"IsoSurfaceValue",
          this->DefaultParameterValues[8] = 0.0; //"NumberOfIterations Diffusion",
          this->DefaultParameterValues[9] = defaultTimeStep ;//"TimeStep",
          this->DefaultParameterValues[10] = 0.0; //"ConductanceParameter",
          this->DefaultParameterValues[11] = 0.0; //"LowerThreshold",
          this->DefaultParameterValues[12] = max; //"UpperThreshold",
          this->DefaultParameterValues[13] = 0.0; // "OutsideValue",
          this->DefaultParameterValues[14] = max; //"InsideValue",
          }
        }
        break;

      case VTK_LEVEL_SET_FAST_MARCHING:
        {
        if (this->DefaultParameterValues.Size() != numberOfParameters )
          {
          this->DefaultParameterValues.SetSize(numberOfParameters);
          }
        this->DefaultParameterValues[0] =  defaultTimeStep;//"TimeStep"  Diffusion
        this->DefaultParameterValues[1] =  0.0; //"NumberOfIterationsOfDiffusion",
        this->DefaultParameterValues[2] =  0.0; //  "ConductanceParameter",
        this->DefaultParameterValues[3] =  0.0; // "Sigma", Gradient filter
        this->DefaultParameterValues[4] =  0.0; // "OutputMinimum" sigmoid
        this->DefaultParameterValues[5] =  1.0; //"OutputMaximum",
        this->DefaultParameterValues[6] =  1.0; // "Alpha",
        this->DefaultParameterValues[7] =  0.0; // "Beta",
        this->DefaultParameterValues[8] =  0.0; // "StoppingValue",FastMarching
        this->DefaultParameterValues[9] =  max; // "UpperThreshold", Threshold
        this->DefaultParameterValues[10] = min; // "LowerThreshold",
        this->DefaultParameterValues[11] = 0.0; //"OutsideValue",
        this->DefaultParameterValues[12] = max; // "InsideValue",
        }
        break;

      case VTK_WATERSHED:
        {
        if (this->DefaultParameterValues.Size() != numberOfParameters)
          {
          this->DefaultParameterValues.SetSize(numberOfParameters);
          }
        this->DefaultParameterValues[0] =  0.0;  //"Level",
        this->DefaultParameterValues[1] =  0.0;   //"Thredhold",
        this->DefaultParameterValues[2] =  defaultTimeStep; //"TimeStep",
        this->DefaultParameterValues[3] =  0.0;   //"ConductanceParameter",
        this->DefaultParameterValues[4] =  0.0;   //"NumberOfIterations",
        }
        break;

      case VTK_KNOWLEDGEBASED:
        {
        if (this->DefaultParameterValues.Size() != numberOfParameters)
          {
          this->DefaultParameterValues.SetSize(numberOfParameters);
          }
        this->DefaultParameterValues[0] =  0.0; // Conductance
        this->DefaultParameterValues[1] =  0.0625;// timestep
        this->DefaultParameterValues[2] =  0.0; //number of iterations
        this->DefaultParameterValues[3] =  min; //min output
        this->DefaultParameterValues[4] =  max; // max output
        this->DefaultParameterValues[5] =  MAX_NUMBER_OF_CLASSES; // init NumberOf Classes
        }
        break;
      }

    return 1;
}

//--------------------------------------------------------------------------
double vtkITKImageSegmentation3D::GetDefaultParameterValue(int i)
{
  if (this->Algorithm < 0 || this->Algorithm >= VTK_NUMBER_OF_ALGORITHMS)
    {
    return INVALID_PARAMETER_VALUE;
   }
  else if (i < 0 ||
           static_cast<size_t>(i) >= this->DefaultParameterValues.Size())
    {
    return INVALID_PARAMETER_VALUE;
    }
  else
    {
    return this->DefaultParameterValues[i];
    }
}

//---------------------------------------------------------------------------
// returns the number of parameters reqired by the algorithm
//---------------------------------------------------------------------------
int vtkITKImageSegmentation3D::GetNumberOfParameters(int algorithm)
{
  size_t i, numberOfParameters = -1;
  const char ** tempNames = NULL;

  if (algorithm < 0 || algorithm > VTK_NUMBER_OF_ALGORITHMS)
    {
    return -1;
    }

  switch (algorithm)
    {
    case VTK_LEVEL_SET_CANNY_EDGE:
      tempNames = vtkCannySegmentationLevelSetParameterName;
      break;
    case VTK_LEVEL_SET_FAST_MARCHING:
      tempNames = vtkFastMarchingSegmentationLevelSetParameterName;
      break;
    case VTK_WATERSHED:
      tempNames = vtkWatershedSegmentationParameterName;
      break;
    case VTK_KNOWLEDGEBASED:
      tempNames =  vtkKnowledgeBasedSegmentationParameterName;
      break;
    }

  if (tempNames)
    {
    numberOfParameters = 0;
    for (i = 0; tempNames[i] != 0; i++)
      {
      numberOfParameters++;
      }
    }

  return numberOfParameters;
}

//---------------------------------------------------------------------------
// returns the names of the Parameters required by the Algorithm
//----------------------------------------------------------------------------
const char* vtkITKImageSegmentation3D::GetParameterName(int algorithm,
                                                        int parameter)
{
  size_t i, numberOfParameters;
  const char **tempNames = NULL;

  if (algorithm < 0 || algorithm > VTK_NUMBER_OF_ALGORITHMS)
    {
    return NULL;
    }

  switch (algorithm)
    {
    case VTK_LEVEL_SET_CANNY_EDGE:
      tempNames = vtkCannySegmentationLevelSetParameterName;
      break;
    case VTK_LEVEL_SET_FAST_MARCHING:
      tempNames =  vtkFastMarchingSegmentationLevelSetParameterName;
      break;
    case VTK_WATERSHED:
      tempNames =  vtkWatershedSegmentationParameterName;
      break;
    case VTK_KNOWLEDGEBASED:
      tempNames =  vtkKnowledgeBasedSegmentationParameterName;
      break;
    }

  numberOfParameters = 0;
  for (i = 0; tempNames[i] != 0; i++)
    {
    numberOfParameters++;
    }

  if (parameter < 0 || static_cast<size_t>(parameter) > numberOfParameters)
    {
    return NULL;
    }

  return tempNames[parameter];
}

//----------------------------------------------------------------------------
//  sets the parameters of the algorithm in VTK
//----------------------------------------------------------------------------
void vtkITKImageSegmentation3D::SetParameter(const char * name, double value)
{
  int i, numberOfParameters;
  const char ** tempNames = NULL;
  this->ParameterChanged = 1;

  if (this->Algorithm < 0 || this->Algorithm >= VTK_NUMBER_OF_ALGORITHMS)
    {
    vtkErrorMacro(<<"SetParameter: set Algorithm before setting parameters");
    return;
    }

  switch (this->Algorithm)
    {
    case VTK_LEVEL_SET_CANNY_EDGE:
      tempNames = vtkCannySegmentationLevelSetParameterName;
      break;
    case VTK_LEVEL_SET_FAST_MARCHING:
      tempNames =  vtkFastMarchingSegmentationLevelSetParameterName;
      break;
    case VTK_WATERSHED:
      tempNames = vtkWatershedSegmentationParameterName;
      break;
    case VTK_KNOWLEDGEBASED:
      tempNames =  vtkKnowledgeBasedSegmentationParameterName;
      break;
    }

  numberOfParameters = 0;
  for (i = 0; tempNames[i] != 0; i++)
    {
    numberOfParameters++;
    }

  if (this->Parameters.Size() !=
      static_cast<unsigned int>(numberOfParameters))
    {
    this->Parameters.SetSize(numberOfParameters);
    if (numberOfParameters > 0)
      {
      this->Parameters.Fill(INVALID_PARAMETER_VALUE);
      }
    }

  for (i = 0; i < numberOfParameters; i++)
    {
    if (strcmp(name, tempNames[i]) == 0)
      {
      if (this->Parameters[i] != value)
        {
        this->Parameters[i] = value;
        }
      this->Modified();
      }
    }
}

//--------------------------------------------------------------------------
// a macro to evaluate an expression for all scalar types
//--------------------------------------------------------------------------

#define vtkTypeCaseMacro(expression) \
      case VTK_FLOAT: { typedef float VTK_TT; expression; } \
        break;

/* The code below if for all types, but we only want float

#define vtkTypeCaseMacro(expression) \
      case VTK_UNSIGNED_LONG: { typedef unsigned long VTK_TT; expression; } \
        break; \
     case VTK_UNSIGNED_SHORT: { typedef unsigned short VTK_TT; expression; } \
        break;\
      case VTK_DOUBLE: { typedef double VTK_TT; expression; } \
        break; \
      case VTK_FLOAT: { typedef float VTK_TT; expression; } \
        break; \
      case VTK_LONG: { typedef long VTK_TT; expression; } \
        break; \
      case VTK_UNSIGNED_LONG: { typedef unsigned long VTK_TT; expression; } \
        break; \
      case VTK_INT: { typedef int VTK_TT; expression; } \
        break; \
      case VTK_UNSIGNED_INT: { typedef unsigned int VTK_TT; expression; } \
        break; \
      case VTK_SHORT: { typedef short VTK_TT; expression; } \
        break; \
      case VTK_UNSIGNED_SHORT: { typedef unsigned short VTK_TT; expression; } \
        break; \
      case VTK_CHAR: { typedef char VTK_TT; expression; } \
        break; \

*/


//---------------------------------------------------------------
//Template function to get the Filters associated
//with the specific Algorithm
//---------------------------------------------------------------
template<class T>
void  vtkITKGetFilters(std::vector<itk::Object*> &Filters,
                       int algorithm,
                       T *)
{
  typedef itk::Image<T, DIMENSIONS> ImageType;
  typedef itk::Image<unsigned short, DIMENSIONS> VtkOutputImageType;
  typedef unsigned short VtkOutputImagePixelType;

  if (Filters.size() > 0)
    {
    Filters.clear();
    }

  switch (algorithm)
    {
    case VTK_LEVEL_SET_CANNY_EDGE:
      {
      typedef float CannyOutputImagePixelType;
      typedef itk::Image< float , DIMENSIONS> CannyOutputImageType;
      typedef itk::GradientAnisotropicDiffusionImageFilter< ImageType,
        ImageType > DiffusionFilterType;
      typedef itk::CannySegmentationLevelSetImageFilter< ImageType, ImageType,
        CannyOutputImagePixelType > CannySegmentationLevelSetImageFilterType;
      typedef itk::BinaryThresholdImageFilter< CannyOutputImageType, ImageType>
        BinaryThresholdingFilterType;

      typename CannySegmentationLevelSetImageFilterType::Pointer
        SmartPointerLevelSet = CannySegmentationLevelSetImageFilterType::New();
      typename DiffusionFilterType::Pointer
        SmartPointerDiffusion = DiffusionFilterType::New();
      typename BinaryThresholdingFilterType::Pointer
        SmartPointerThresholder = BinaryThresholdingFilterType::New();

      SmartPointerLevelSet->Register();
      SmartPointerDiffusion->Register();
      SmartPointerThresholder->Register();

      Filters.push_back((itk::Object*) SmartPointerLevelSet.GetPointer());
      Filters.push_back((itk::Object*) SmartPointerDiffusion.GetPointer());
      Filters.push_back((itk::Object*) SmartPointerThresholder.GetPointer());
      }
      break;

    case VTK_LEVEL_SET_FAST_MARCHING:
      {
      //Filter types
      typedef itk::GradientAnisotropicDiffusionImageFilter
        < ImageType, ImageType > DiffusionFilterType;
      typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< ImageType,
        ImageType> GradientMagnitudeFilterType;
      typedef itk::SigmoidImageFilter< ImageType , ImageType > SigmoidFilterType;
      typedef itk::FastMarchingImageFilter<
        ImageType, ImageType > FastMarchingFilterType;
      typedef itk::BinaryThresholdImageFilter<
        ImageType, ImageType > ThresholdFilterType;
      //Setting Filters
      typename DiffusionFilterType:: Pointer
        diffusionSmartPointer = DiffusionFilterType::New();
      typename GradientMagnitudeFilterType:: Pointer
        gradientSmartPointer = GradientMagnitudeFilterType::New();
      typename SigmoidFilterType::Pointer
        sigmoidSmartPointer = SigmoidFilterType::New();
      typename FastMarchingFilterType::Pointer
        fastmarchingSmartPointer =  FastMarchingFilterType::New();
      typename ThresholdFilterType::Pointer
        thresholdSmartPointer = ThresholdFilterType::New();

      //Increasing the reference count
      diffusionSmartPointer->Register();
      gradientSmartPointer->Register();
      sigmoidSmartPointer->Register();
      fastmarchingSmartPointer->Register();
      thresholdSmartPointer->Register();

      Filters.push_back(dynamic_cast<itk::Object*>(
                          diffusionSmartPointer.GetPointer()));
      Filters.push_back(dynamic_cast<itk::Object*>(
                          gradientSmartPointer.GetPointer()));
      Filters.push_back(dynamic_cast<itk::Object*>(
                          sigmoidSmartPointer.GetPointer()));
      Filters.push_back(dynamic_cast<itk::Object*>(
                          fastmarchingSmartPointer.GetPointer()));
      Filters.push_back(dynamic_cast<itk::Object*>(
                          thresholdSmartPointer.GetPointer()));
      }
      break;

    case VTK_WATERSHED:
      {
      typedef itk::GradientAnisotropicDiffusionImageFilter
        < ImageType, ImageType > DiffusionFilterType;
      typedef itk::WatershedImageFilter<ImageType> WatershedFilterType;
      typedef itk::GradientMagnitudeImageFilter< ImageType,
        ImageType> GradientMagnitudeFilterType;

      typename DiffusionFilterType::Pointer
        SmartPointerDiffusion = DiffusionFilterType::New();
      typename GradientMagnitudeFilterType::Pointer
        SmartPointerGradientMagnitude = GradientMagnitudeFilterType::New();
      typename WatershedFilterType::Pointer
        SmartPointerWatershed =  WatershedFilterType::New();

      SmartPointerDiffusion->Register();
      SmartPointerGradientMagnitude->Register();
      SmartPointerWatershed->Register();

      Filters.push_back((itk::Object*)SmartPointerWatershed.GetPointer());
      Filters.push_back((itk::Object*)SmartPointerDiffusion.GetPointer());
      Filters.push_back((itk::Object*)SmartPointerGradientMagnitude.GetPointer());
      }
      break;

    case VTK_KNOWLEDGEBASED:
      {
      typedef itk::VectorImage<T, DIMENSIONS> VectorImageType;
      typedef short LabelPixelType;
      typedef itk::Image<LabelPixelType, DIMENSIONS> LabelImageType;

      typedef itk::BayesianClassifierInitializationImageFilter<
        ImageType> MembershipImageGeneratorType;
      typedef itk::BayesianClassifierImageFilter<
        VectorImageType,LabelPixelType > ClassifierFilterType;
      typedef itk::Image<double, DIMENSIONS> PosteriorImageType;
      typedef itk::GradientAnisotropicDiffusionImageFilter<
       PosteriorImageType, PosteriorImageType >  SmoothingFilterType;
      typedef itk::RescaleIntensityImageFilter< LabelImageType,
        VtkOutputImageType >    RescaleFilterType;

      typename MembershipImageGeneratorType::Pointer smartPointerMembershipImageGenerator
        = MembershipImageGeneratorType::New();
      typename ClassifierFilterType::Pointer smartPointerClassifier = ClassifierFilterType::New();
      typename SmoothingFilterType::Pointer smartPointerSmoother = SmoothingFilterType::New();
      typename RescaleFilterType::Pointer smartPointerRescale = RescaleFilterType::New();

      smartPointerMembershipImageGenerator->Register();
      smartPointerClassifier->Register();
      smartPointerSmoother->Register();
      smartPointerRescale->Register();

      Filters.push_back((itk::Object*)smartPointerMembershipImageGenerator.GetPointer());
      Filters.push_back((itk::Object*)smartPointerClassifier.GetPointer());
      Filters.push_back((itk::Object*)smartPointerSmoother.GetPointer());
      Filters.push_back((itk::Object*)smartPointerRescale.GetPointer());
      }
      break;
    }// end of switch
}

//---------------------------------------------------------------
//Template function to get the set of Filters according
//to the Type of the Algorithm
//---------------------------------------------------------------
template<class T>
void  vtkITKSetFilterParameters(std::vector< itk::Object* > Filters,
                                itk::Array< double > Parameters,
                                itk::Array< double > DefaultParameterValues,
                                int algorithm,
                                T *)
{
  typedef itk::Image<T, DIMENSIONS> ImageType;
  typedef itk::Image<unsigned short, DIMENSIONS> VtkOutputImageType;
  typedef unsigned short VtkOutputImagePixelType;

  switch (algorithm)
    {
    case VTK_LEVEL_SET_FAST_MARCHING:
      {
      typedef itk::GradientAnisotropicDiffusionImageFilter
        < ImageType, ImageType > DiffusionFilterType;
      typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< ImageType,
        ImageType> GradientMagnitudeFilterType;
      typedef itk::SigmoidImageFilter< ImageType, ImageType > SigmoidFilterType;
      typedef itk::FastMarchingImageFilter<
        ImageType, ImageType > FastMarchingFilterType;
      typedef itk::BinaryThresholdImageFilter<
        ImageType, ImageType > ThresholdFilterType;
      typedef typename ThresholdFilterType::OutputPixelType
        ThresholdOutputPixelType;
      typedef typename ThresholdFilterType::InputPixelType
        ThresholdInputPixelType;

      DiffusionFilterType* diffusion
        = dynamic_cast<DiffusionFilterType*>( Filters.at(0));
      GradientMagnitudeFilterType * gradient
        = dynamic_cast<GradientMagnitudeFilterType *>(Filters.at(1));
      SigmoidFilterType* sigmoid
        = dynamic_cast<SigmoidFilterType*>(Filters.at(2));
      FastMarchingFilterType * fastmarching
        = dynamic_cast<FastMarchingFilterType*>(Filters.at(3));
      ThresholdFilterType * threshold
        = dynamic_cast<ThresholdFilterType*>(Filters.at(4));

      for ( int i = 0; i<= 12; i ++ )
        {
        if ( Parameters[i] == INVALID_PARAMETER_VALUE)
          {
          Parameters[i] =  DefaultParameterValues[i];
          }
        }
      diffusion->SetTimeStep( Parameters[0] );
      diffusion->SetNumberOfIterations(
        static_cast<unsigned int>(Parameters[1]) );
      diffusion->SetConductanceParameter( Parameters[2] );
      gradient->SetSigma( Parameters[3] );
      sigmoid->SetOutputMinimum( static_cast<unsigned short>(Parameters[4]) );
      sigmoid->SetOutputMaximum(static_cast<unsigned short>(Parameters[5]) );
      sigmoid->SetAlpha( Parameters[6] );
      sigmoid->SetBeta( Parameters[7] );
      fastmarching->SetStoppingValue(static_cast< int >(Parameters[8]));
      threshold->SetUpperThreshold(
        static_cast<ThresholdInputPixelType>(Parameters[9]));
      threshold->SetLowerThreshold(
        static_cast<ThresholdInputPixelType>(Parameters[10]));
      threshold->SetOutsideValue(
        static_cast<ThresholdOutputPixelType>(Parameters[11]));
      threshold->SetInsideValue(
        static_cast<ThresholdOutputPixelType>(Parameters[12]));
      }
      break;

    case VTK_LEVEL_SET_CANNY_EDGE:
      {
      for (int i = 0; i <= 14; i++)
        {
        if (Parameters[i] == INVALID_PARAMETER_VALUE)
          {
          Parameters[i] =  DefaultParameterValues[i];
          }
        }
      typedef float CannyOutputImagePixelType;
      typedef itk::Image< float , DIMENSIONS> CannyOutputImageType;
      typedef itk::GradientAnisotropicDiffusionImageFilter< ImageType,
        ImageType > DiffusionFilterType;
      typedef itk::CannySegmentationLevelSetImageFilter< ImageType, ImageType,
        CannyOutputImagePixelType> CannySegmentationLevelSetImageFilterType;
      typedef itk::BinaryThresholdImageFilter< CannyOutputImageType, ImageType >
        BinaryThresholdFilterType;
      typedef typename BinaryThresholdFilterType::OutputPixelType
        ThresholdOutputPixelType;
      typedef typename BinaryThresholdFilterType::InputPixelType
        ThresholdInputPixelType;

      CannySegmentationLevelSetImageFilterType* cannyPointer =
        dynamic_cast<CannySegmentationLevelSetImageFilterType*>(Filters.at(0));
      DiffusionFilterType* diffusionPointer =
        dynamic_cast<DiffusionFilterType*>(Filters.at(1));
      BinaryThresholdFilterType *thresholdPointer =
        dynamic_cast<BinaryThresholdFilterType *>(Filters.at(2));

      cannyPointer->SetThreshold(Parameters[0]);
      cannyPointer->SetVariance(Parameters[1]);
      cannyPointer->SetCurvatureScaling(Parameters[2]);
      cannyPointer->SetPropagationScaling(Parameters[3]);
      cannyPointer->SetAdvectionScaling(Parameters[4]);
      cannyPointer->SetMaximumRMSError(Parameters[5]);
      cannyPointer->SetNumberOfIterations((int)Parameters[6]);
      cannyPointer->SetIsoSurfaceValue(Parameters[7]);
      diffusionPointer->SetNumberOfIterations(
        static_cast<unsigned int> (Parameters[8]));
      diffusionPointer->SetTimeStep(Parameters[9]);
      diffusionPointer->SetConductanceParameter(Parameters[10]);
      thresholdPointer->SetLowerThreshold(static_cast<float>(Parameters[11]));
      thresholdPointer->SetUpperThreshold(
        static_cast<ThresholdInputPixelType>(Parameters[12]));
      thresholdPointer->SetOutsideValue(
        static_cast<ThresholdOutputPixelType>(Parameters[13]));
      thresholdPointer->SetInsideValue(
        static_cast<ThresholdOutputPixelType>(Parameters[14]));
      }
      break;

    case VTK_WATERSHED:
      {
      for (int i = 0; i <= 4; i++)
        {
        if (Parameters[i] == INVALID_PARAMETER_VALUE)
          {
          Parameters[i] =  DefaultParameterValues[i];
          }
        }
      typedef itk::GradientAnisotropicDiffusionImageFilter
        < ImageType, ImageType > DiffusionFilterType;
      typedef itk::GradientMagnitudeImageFilter< ImageType ,
        ImageType> GradientMagnitudeFilterType;
      typedef itk::WatershedImageFilter< ImageType > WatershedFilterType;

      WatershedFilterType* watershedPointer =
        dynamic_cast<WatershedFilterType*>(Filters.at(0));
      DiffusionFilterType* diffusionPointer =
        dynamic_cast<DiffusionFilterType*>(Filters.at(1));

      watershedPointer->SetLevel(Parameters[0]);
      watershedPointer->SetThreshold(Parameters[1]);
      diffusionPointer->SetTimeStep(Parameters[2]);
      diffusionPointer->SetConductanceParameter(Parameters[3]);
      diffusionPointer->SetNumberOfIterations(static_cast<int>(Parameters[4]));
      }
      break;

    case VTK_KNOWLEDGEBASED:
      {
      int n = Parameters.Size();
      for (int i = 0; i <= n; i++)
        {
        if (Parameters[i] == INVALID_PARAMETER_VALUE)
          {
          Parameters[i] =  DefaultParameterValues[i];
          }
        }

      typedef itk::VectorImage<T, DIMENSIONS> VectorImageType;
      typedef short LabelPixelType;
      typedef itk::Image<LabelPixelType, DIMENSIONS> LabelImageType;

      typedef itk::BayesianClassifierInitializationImageFilter<
        ImageType> MembershipImageGeneratorType;
      typedef itk::BayesianClassifierImageFilter<
        VectorImageType,LabelPixelType > ClassifierFilterType;
      typedef itk::Image<double, DIMENSIONS> PosteriorImageType;
      typedef itk::GradientAnisotropicDiffusionImageFilter<
       PosteriorImageType, PosteriorImageType >  SmoothingFilterType;
      typedef itk::RescaleIntensityImageFilter< LabelImageType,
        VtkOutputImageType >    RescaleFilterType;

      MembershipImageGeneratorType* membershipImageGenerator
        = dynamic_cast<MembershipImageGeneratorType*>(Filters.at(0));
      ClassifierFilterType* classifier = dynamic_cast<ClassifierFilterType*>(Filters.at(1));
      SmoothingFilterType* smoother = dynamic_cast<SmoothingFilterType*>(Filters.at(2));
      RescaleFilterType* rescale = dynamic_cast<RescaleFilterType*>(Filters.at(3));

      membershipImageGenerator->SetNumberOfClasses(static_cast<int>(Parameters[5]));
      smoother->SetTimeStep(Parameters[0]);
      smoother->SetNumberOfIterations(static_cast<int>(Parameters[1]));
      smoother->SetConductanceParameter(Parameters[2]);

      classifier->SetNumberOfSmoothingIterations(static_cast<int>(Parameters[1]));

      rescale->SetOutputMinimum(
        static_cast<VtkOutputImagePixelType>(Parameters[3]));
      rescale->SetOutputMaximum(
        static_cast<VtkOutputImagePixelType>(Parameters[4]));
      }
      break;
    }
}

//--------------------------------------------------------------------------
// Template function to Export image from VTK to ITK
//--------------------------------------------------------------------------
template<class T> static void
vtkITK3DConnectVTKITKPipelines(vtkImageExport *vtkImageExporter,
                               vtkImageExport *vtkLabelImageExporter,
                               itk::Object    *&itkImageImporter,
                               itk::Object    *&itkLabelImageImporter,
                               T *)
{
  typedef itk::Image<T, DIMENSIONS> ImageType;

  typename itk::VTKImageImport< ImageType >::Pointer imageSmartPointer
    = itk::VTKImageImport< ImageType >::New();
  imageSmartPointer->Register();
  typename itk::VTKImageImport< ImageType >* imagePointer
    = imageSmartPointer.GetPointer();
  itkImageImporter = imagePointer;

  ::ConnectVTKToITK(vtkImageExporter, imagePointer);

  if(vtkLabelImageExporter != NULL)
    {
    typename itk::VTKImageImport< ImageType >::Pointer labelSmartPointer
      = itk::VTKImageImport< ImageType >::New();
    labelSmartPointer->Register();
    typename itk::VTKImageImport< ImageType >* labelPointer
      = labelSmartPointer.GetPointer();
    itkLabelImageImporter = labelPointer;

    ::ConnectVTKToITK(vtkLabelImageExporter, labelPointer);
    }
}

//--------------------------------------------------------------------------
// Template function to export Image from ITK to VTK
//--------------------------------------------------------------------------
template<class T>
static void vtkITK3DConnectITKVTKPipelines(vtkImageImport *vtkImageImporter,
                                           itk::Object    *&itkImageExporter,
                                           T *)
{
  typedef itk::Image< unsigned short , DIMENSIONS > VtkOutputImageType;
  typename itk::VTKImageExport<VtkOutputImageType >::Pointer
    segmentSmartPointer = itk::VTKImageExport< VtkOutputImageType >::New();
  segmentSmartPointer->Register();

  typename itk::VTKImageExport< VtkOutputImageType >* segmentPointer
    =  segmentSmartPointer.GetPointer();
  itkImageExporter = segmentPointer;

  ::ConnectITKToVTK(segmentPointer, vtkImageImporter);
}

//--------------------------------------------------------------------------
// Template function to import seed from ITK to VTK
//--------------------------------------------------------------------------
template<class T>
static void vtkITKSetSeeds(vtkPoints *vtkseeds,
                    int algorithm,
                    std::vector<itk::Object*>Filters,
                    T *)
{
  if (algorithm == VTK_LEVEL_SET_FAST_MARCHING)
    {
    typedef itk::Image< T , DIMENSIONS> ImageType;
    typedef itk::FastMarchingImageFilter<
      ImageType, ImageType > FastMarchingFilterType;
    typename ImageType::IndexType itkSeed;
    typedef typename FastMarchingFilterType::NodeContainer  NodeContainer;
    typedef typename FastMarchingFilterType::NodeType   NodeType;
    typename NodeContainer::Pointer seeds = NodeContainer::New();
    int numberofseeds =0;

    if (vtkseeds != NULL)
      {
      numberofseeds =  vtkseeds->GetNumberOfPoints();
      }
    else
      {
      cerr<<"No seed has been entered\n";
      return;
      }

    double temp[DIMENSIONS];
    const unsigned short seedValue = 0;
    NodeType node;
    FastMarchingFilterType * fastMarching
      = dynamic_cast<FastMarchingFilterType*>(Filters.at(3));
    int i,j;

    if (numberofseeds > 0)
      {
      seeds->Initialize();
      }

    for (i = 0; i < numberofseeds; i++)
      {
      vtkseeds->GetPoint(i, temp);
      for (j = 0; j < DIMENSIONS; j++)
        {
        itkSeed[j] = static_cast<long >(temp[j]);
        }
      node.SetValue(seedValue);
      node.SetIndex(itkSeed);
      seeds->InsertElement(i , node);
      }

    if ( numberofseeds > 0 )
      {
      fastMarching->SetTrialPoints(seeds);
      }
    }
}

//-------------------------------------------------------------------
// This Template Function Performs the Segmentation Steps
//-------------------------------------------------------------------
template< class T >
static void vtkITKSegment(std::vector<itk::Object*>Filters,
                   int algorithm,
                   itk::Object* itkImageImporter,
                   itk::Object* itkLabelImageImporter,
                   itk::Object* itkImageExporter,
                   T *)
{
  typedef itk::Image< T , DIMENSIONS> ImageType;
  typedef itk::Image< unsigned short , DIMENSIONS > VtkOutputImageType;
  typedef itk::VTKImageExport< VtkOutputImageType > VTKImageExportType;

  switch(algorithm)
    {
    case VTK_LEVEL_SET_CANNY_EDGE:
      {
      typedef float CannyOutputImagePixelType;
      typedef itk::Image<float, DIMENSIONS > CannyOutputImageType;
      typedef itk::GradientAnisotropicDiffusionImageFilter<
        ImageType, ImageType > DiffusionFilterType;
      typedef itk::CannySegmentationLevelSetImageFilter<ImageType, ImageType,
        CannyOutputImagePixelType > CannySegmentationLevelSetImageFilterType;
      typedef itk::BinaryThresholdImageFilter< CannyOutputImageType,
        ImageType > BinaryThresholdFilterType;

      CannySegmentationLevelSetImageFilterType* canny =
        (dynamic_cast<CannySegmentationLevelSetImageFilterType *>
         (Filters.at(0)));
      DiffusionFilterType *diffusion =
        (dynamic_cast<DiffusionFilterType *> (Filters.at(1)));
      BinaryThresholdFilterType * threshold =
        (dynamic_cast<BinaryThresholdFilterType *>(Filters.at(2)));
      typedef  itk::CastImageFilter< ImageType, VtkOutputImageType > CastType;
      typename CastType::Pointer caster = CastType::New();

      diffusion->SetInput((dynamic_cast<itk::VTKImageImport< ImageType >*>
                           (itkImageImporter))->GetOutput());
      canny->SetInitialImage((dynamic_cast<itk::VTKImageImport< ImageType >*>
                              (itkLabelImageImporter))->GetOutput());
      canny->SetFeatureImage(diffusion->GetOutput());
      threshold->SetInput(canny->GetOutput());
      caster->SetInput(threshold->GetOutput());
      caster->Update();
      (dynamic_cast<VTKImageExportType *>(itkImageExporter))
        ->SetInput(caster->GetOutput());
      }
      break;

    case VTK_LEVEL_SET_FAST_MARCHING:
      {
      typedef itk::GradientAnisotropicDiffusionImageFilter
        < ImageType, ImageType > DiffusionFilterType;
      typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<
        ImageType, ImageType> GradientMagnitudeFilterType;
      typedef itk::SigmoidImageFilter<
        ImageType , ImageType > SigmoidFilterType;
      typedef itk::FastMarchingImageFilter<
        ImageType, ImageType > FastMarchingFilterType;
      typedef itk::BinaryThresholdImageFilter<
        ImageType, ImageType > ThresholdFilterType;
      typedef typename ThresholdFilterType::OutputPixelType
        ThresholdOutputPixelType;
      typedef typename ThresholdFilterType::InputPixelType
        ThresholdInputPixelType;
      typedef typename itk::VTKImageImport< ImageType > VTKImageImportType;

      DiffusionFilterType* diffusion
        = dynamic_cast<DiffusionFilterType*>( Filters.at(0));
      GradientMagnitudeFilterType * gradient
          = dynamic_cast<GradientMagnitudeFilterType *>(Filters.at(1));
      SigmoidFilterType* sigmoid
        = dynamic_cast<SigmoidFilterType*>(Filters.at(2));
      FastMarchingFilterType * fastmarching
        = dynamic_cast<FastMarchingFilterType*>(Filters.at(3));
      ThresholdFilterType * threshold
        = dynamic_cast<ThresholdFilterType*>(Filters.at(4));

      diffusion->SetInput((dynamic_cast<VTKImageImportType*>
                           (itkImageImporter))->GetOutput());
      gradient->SetInput(diffusion->GetOutput());
      sigmoid->SetInput(gradient->GetOutput());
      fastmarching->SetInput(sigmoid->GetOutput());
      threshold->SetInput(fastmarching->GetOutput());

      typedef  itk::CastImageFilter< ImageType, VtkOutputImageType > CastType;
      typename CastType::Pointer caster = CastType::New();
      caster->SetInput(threshold->GetOutput());
      caster->Update();
      (dynamic_cast<VTKImageExportType *>(itkImageExporter))
        ->SetInput(caster->GetOutput());
      }
      break;

    case VTK_WATERSHED:
      {
      typedef itk::GradientAnisotropicDiffusionImageFilter
        < ImageType, ImageType > DiffusionFilterType;
      typedef itk::GradientMagnitudeImageFilter< ImageType,
          ImageType> GradientMagnitudeFilterType;
      typedef itk::Image< T  ,DIMENSIONS > OutputImageType;
      typedef itk::WatershedImageFilter<ImageType> WatershedFilterType;
      typedef unsigned long WatershedOutputPixelType;
      typedef itk::Image< WatershedOutputPixelType ,
        DIMENSIONS > WatershedOutputImageType;
      typedef  itk::CastImageFilter< WatershedOutputImageType,
        VtkOutputImageType > CastType;

      GradientMagnitudeFilterType* gradient =
        dynamic_cast<GradientMagnitudeFilterType *>(Filters.at(2));
      DiffusionFilterType*  diffusion =
        dynamic_cast<DiffusionFilterType *>(Filters.at(1));
      WatershedFilterType* watershed =
        dynamic_cast<WatershedFilterType *>(Filters.at(0));

      diffusion->SetInput((dynamic_cast<itk::VTKImageImport< ImageType >*>
                           (itkImageImporter))->GetOutput());
      gradient->SetInput(  diffusion->GetOutput() );
      watershed->SetInput( gradient->GetOutput()  );

      typename CastType::Pointer caster = CastType::New();
      caster->SetInput(watershed->GetOutput());
      caster->Update();
      (dynamic_cast<VTKImageExportType *>(itkImageExporter))
        ->SetInput(caster->GetOutput());
      }
      break;

    case VTK_KNOWLEDGEBASED:
      {
      typedef itk::VectorImage<T, DIMENSIONS> VectorImageType;
      typedef short LabelPixelType;
      typedef itk::Image<LabelPixelType, DIMENSIONS> LabelImageType;

      typedef itk::BayesianClassifierInitializationImageFilter<
        ImageType> MembershipImageGeneratorType;
      typedef itk::BayesianClassifierImageFilter<
        VectorImageType,LabelPixelType > ClassifierFilterType;
      typedef itk::Image<double, DIMENSIONS> PosteriorImageType;
      typedef itk::GradientAnisotropicDiffusionImageFilter<
        PosteriorImageType, PosteriorImageType >  SmoothingFilterType;
      typedef itk::RescaleIntensityImageFilter< LabelImageType,
        VtkOutputImageType >    RescaleFilterType;

      MembershipImageGeneratorType* membershipImageGenerator
        = dynamic_cast<MembershipImageGeneratorType*>(Filters.at(0));
      ClassifierFilterType* classifier = dynamic_cast<ClassifierFilterType*>(Filters.at(1));
      SmoothingFilterType* smoother = dynamic_cast<SmoothingFilterType*>(Filters.at(2));
      RescaleFilterType* rescale = dynamic_cast<RescaleFilterType*>(Filters.at(3));

      membershipImageGenerator->SetInput(
             (dynamic_cast<itk::VTKImageImport< ImageType >*>
              (itkImageImporter))->GetOutput());

      classifier->SetInput(membershipImageGenerator->GetOutput());
      classifier->SetSmoothingFilter(smoother);
      //classifier->GetOutput()->Update();

      rescale->SetInput(classifier->GetOutput());
      rescale->Update();
      (dynamic_cast<VTKImageExportType *>(itkImageExporter))
        ->SetInput(rescale->GetOutput());
      }
      break;
    }
}

//---------------------------------------------------------------------------
// The main module to set up the segmentation framework and connect
// pieces together and do Segmentation
int vtkITKImageSegmentation3D::SetVTKITKPipelineConnection()
{
  void *dummy = NULL;
  vtkPoints *points = 0;
  vtkPoints *voxelPoints = 0;

  if (this->GetMTime() > this->ExecuteTime)
    {
    int ret = this->Initialize();
    if (ret == 1)
      {
      vtkErrorMacro(<<"vtkITKImageSegmentation3D::"
                    <<"SetVTKITKPipelineConnection - Initialize failed\n");
      return 0;
      }

    this->InvokeEvent(vtkCommand::StartEvent,NULL);
    // reset Abort flag
    this->AbortExecute = 0;
    this->Progress = 0.0;

    /////////////////////////////////////////////////////////////
    // Get the Filters to be used
    switch (this->FeatureImageDataType)
      {
      vtkTypeCaseMacro( vtkITKGetFilters(
                          this->Filters,
                          this->Algorithm,
                          (VTK_TT *)dummy) );
      }

    ////////////////////////////////////////////////////////////
    // Setup the Parameters of the filters  <ITK>
    switch (this->FeatureImageDataType)
      {
      vtkTypeCaseMacro( vtkITKSetFilterParameters(
                          this->Filters,
                          this->Parameters,
                          this->DefaultParameterValues,
                          this->Algorithm,
                          (VTK_TT *)dummy) );
      }

    ////////////////////////////////////////////////////////////
    //Setup the list of seeds in ITK
    if (this->GetNeedsSeedPoints(this->Algorithm))
      {
      points = this->GetSeedPoints()->GetPoints();
      points->Register(this);
      voxelPoints = vtkPoints::New();
      int n = points->GetNumberOfPoints();
      double w[3];
      double v[3];
      voxelPoints->SetNumberOfPoints(n);
      for(int i= 0; i<n;i++)
        {
        points->GetPoint(i, w);
        this->ConvertWorldToVoxel(w,v);
        voxelPoints->InsertPoint(i, v);
        }

      switch (this->FeatureImageDataType)
        {
        vtkTypeCaseMacro( vtkITKSetSeeds(
                        voxelPoints,
                        this->Algorithm,
                        this->Filters,
                        (VTK_TT *)dummy));
        }
      points->Delete();
      voxelPoints->Delete();
      }

    ///////////////////////////////////////////////////////
    // Connect vtk and itk pipelines

    vtkImageExport* vtkImageExporter = vtkImageExport::New();
    vtkImageExporter->SetInput( this->GetFeatureImage() );
    vtkImageExport* vtkLabelImageExporter = vtkImageExport::New();

    if(this->GetNeedsLabelImage(this->Algorithm))
      {
      vtkLabelImageExporter->SetInput( this->GetLabelImage() );
      }

    switch (this->FeatureImageDataType)
      {
      vtkTypeCaseMacro( vtkITK3DConnectVTKITKPipelines(
                          vtkImageExporter,
                          vtkLabelImageExporter,
                          this->ITKImageImporter,
                          this->ITKLabelImageImporter,
                          (VTK_TT *)dummy) );

      }

    ///////////////////////////////////////////////////////
    // Connect itk and vtk pipelines

    switch (this->FeatureImageDataType)
      {
      vtkTypeCaseMacro( vtkITK3DConnectITKVTKPipelines(
                          this->VTKImageImporter,
                          this->ITKImageExporter,
                          (VTK_TT *)dummy) );

      }

    /////////////////////////////////////////////////////
    // Performs the segmentation

    switch (this->FeatureImageDataType)
      {
      vtkTypeCaseMacro(vtkITKSegment(
                          this->Filters,
                          this->Algorithm,
                          this->ITKImageImporter,
                          this->ITKLabelImageImporter,
                          this->ITKImageExporter,
                          (VTK_TT *)dummy) );
      }

    vtkImageData *importOutput = this->VTKImageImporter->GetOutput();

    importOutput->Update();

    // Create the stencil output
    vtkImageToImageStencil *imageToStencil = vtkImageToImageStencil::New();

    imageToStencil->SetInput(importOutput);
    imageToStencil->ThresholdByUpper(this->StencilThreshold);
    //cout<<"Calling Update on imageToStencil\n";
    // imageToStencil->Update();

    vtkImageStencilData *stencilOutput = imageToStencil->GetOutput();
     // Pass the data to our first output
    vtkInformation *outInfo = this->GetExecutive()->GetOutputInformation(0);
    vtkImageData *output = vtkImageData::SafeDownCast(
      outInfo->Get(vtkDataObject::DATA_OBJECT()));

    output->SetExtent(importOutput->GetExtent());
    output->GetPointData()->PassData(importOutput->GetPointData());

    // Pass the stencil to our second output
    vtkInformation *outInfo2 = this->GetExecutive()->GetOutputInformation(1);
    vtkImageStencilData *output2 = vtkImageStencilData::SafeDownCast(
      outInfo2->Get(vtkDataObject::DATA_OBJECT()));

    output2->SetExtent(stencilOutput->GetExtent());
    output2->InternalImageStencilDataCopy(stencilOutput);
    cout<<" output updated\n";
    this->ExecuteTime.Modified();
    if ( !this->AbortExecute )
      {
      this->UpdateProgress(1.0);
      }
    this->InvokeEvent(vtkCommand::EndEvent,NULL);
    }

  return 1;
}

//----------------------------------------------------------------------------
void vtkITKImageSegmentation3D::ConvertWorldToVoxel(double w[3],
                                                    double v[3])
{
  double imageOrigin[3], imageSpacing[3];
  this->GetInput()->GetOrigin(imageOrigin);
  this->GetInput()->GetSpacing(imageSpacing);
  v[0] = double(vtkMath::Round(((w[0]-imageOrigin[0])/imageSpacing[0])));
  v[1] = double(vtkMath::Round(((w[1]-imageOrigin[1])/imageSpacing[1])));
  v[2] = double(vtkMath::Round(((w[2]-imageOrigin[2])/imageSpacing[2])));
  return;
}

//----------------------------------------------------------------------------
// this is not an initialization routine rather a rountine to check
// for correct data type and components
int vtkITKImageSegmentation3D::Initialize()
{
  vtkImageData* image = this->GetFeatureImage();

  if ( image == NULL )
    {
    vtkErrorMacro(<<"Image is not present\n");
    return 1;
    }

  if (this->GetNeedsLabelImage(this->Algorithm))
    {
    if (this->GetLabelImage() == NULL)
      {
      vtkErrorMacro(<<"Label Image is not present\n");
      return 1;
      }
    }

  if (this->GetNeedsSeedPoints(this->Algorithm))
    {
    if( this->GetSeedPoints() == NULL )
      {
      vtkErrorMacro(<<"Seeds are not set!\n");
      return 1;
      }
    }

  if( this->Algorithm < 0 || this->Algorithm >= VTK_NUMBER_OF_ALGORITHMS)
    {
    vtkErrorMacro(<<"Algorithm is not supported\n");
    return 1;
    }

  // Initialization OK.
  return 0;
}
