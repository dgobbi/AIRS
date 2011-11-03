/*=========================================================================

  Program:   Atamai Classes for VTK
  Module:    $RCSfile: vtkITKKnowledgeBasedClassification.cxx,v $
  Creator:   Piali Das <pdas@atamai.com>
  Language:  C++
  Author:    $Author: dgobbi $
  Date:      $Date: 2007/05/07 22:42:30 $
  Version:   $Revision: 1.3 $

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

#include "vtkITKKnowledgeBasedClassification.h"

// VTK header files

#include "vtkPoints.h"
#include "vtkPointSet.h"
#include "vtkTransform.h"
#include "vtkTimerLog.h"
#include "vtkExecutive.h"
#include "vtkObjectFactory.h"
#include "vtkImageData.h"
#include "vtkImageImport.h"
#include "vtkImageExport.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkCommand.h"

#include "math.h"
#include "itkImage.h"
#include "itkVector.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageConstIterator.h"
#include "itkJoinImageFilter.h"
#include "itkScalarImageKmeansImageFilter.h"

#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageCastVectorIndexSelectionFilter.h"
#include "itkBayesianClassifierImageFilter.h"
#include "itkVectorImage.h"
#include "itkBayesianClassifierInitializationImageFilter.h"

#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkMaximumDecisionRule.h"
#include "itkRescaleIntensityImageFilter.h"



#include "vtkITKImageAlgorithm.h"

// itk image import/export files
#include "itkImage.h"
#include "itkVTKImageImport.h"
#include "itkVTKImageExport.h"
#include "ConnectVTKITK.h"
#include "itkCastImageFilter.h"

#include "itkCommand.h"

// Parameters
const char * KnowldgeBased[] = {
  "InitialMean1",
  "InitialMean2",
  "InitialMean3",
  "Conductance",
  "TimeStep",
  "NumberOfIterations",
  "OutputMinimum",
  "OutputMaximum",
  0
};

//----------------------------------------------------------------------------
// Implements the New method for the class
//----------------------------------------------------------------------------
vtkITKKnowledgeBasedClassification* vtkITKKnowledgeBasedClassification::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret =
    vtkObjectFactory::CreateInstance("vtkITKKnowledgeBasedClassification");

  if(ret)
    {
    return (vtkITKKnowledgeBasedClassification*)ret;
    }

  // If the factory was unable to create the object, then create it here.
  return new vtkITKKnowledgeBasedClassification;
}


//----------------------------------------------------------------------------
//CONSTRUCTOR
//----------------------------------------------------------------------------
vtkITKKnowledgeBasedClassification::vtkITKKnowledgeBasedClassification()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(0);

  this->NumberOfClasses = 2;
  this->Mean = new double[2];
  this->Variance = new double[2];
  for (int i = 0; i <this->NumberOfClasses ; i++)
    {
    this->Mean[i] = 0.0;
    this->Variance[i] = 0.0;
    }

  this->ImageDimension = 3;
  this->ImageDataType = 10;
}

//----------------------------------------------------------------------------
//DESTRUCTTOR
//----------------------------------------------------------------------------
vtkITKKnowledgeBasedClassification::~vtkITKKnowledgeBasedClassification()
{

}

//----------------------------------------------------------------------------
// PRINT SELF
//---------------------------------------------------------------------------
 void vtkITKKnowledgeBasedClassification::PrintSelf(ostream& os, vtkIndent indent)
{
  this->PrintSelf(os,indent);
}


//----------------------------------------------------------------------------
void vtkITKKnowledgeBasedClassification::SetNumberOfClasses(const int n)
{
  if( n>0 && this->NumberOfClasses != n)
    {
    this->NumberOfClasses = n;
    this->Mean = new double[n];
    this->Variance = new double[n];
    for (int i = 0; i < n; i++)
      {
      this->Mean[i] = 0.0;
      this->Variance[i] = 0.0;
      }
    }
}
//----------------------------------------------------------------------------
void vtkITKKnowledgeBasedClassification::SetMean(int index, double m)
{
  if( index >= 0 && index < this->NumberOfClasses )
    {
    this->Mean[index] = m;
    }
}
//----------------------------------------------------------------------------
double vtkITKKnowledgeBasedClassification::GetMean(int index)
{
  if( index >= 0 && index < this->NumberOfClasses )
    {
    return this->Mean[index];
    }

  return 0.0;
}
//----------------------------------------------------------------------------
void vtkITKKnowledgeBasedClassification::SetVariance(int index, double v)
{
  if( index >= 0 && index < this->NumberOfClasses )
    {
    this->Variance[index] = v;
    }
}
//----------------------------------------------------------------------------
double vtkITKKnowledgeBasedClassification::GetVariance(int index)
{
  if( index >= 0 && index < this->NumberOfClasses )
    {
    return this->Variance[index];
    }

  return 0.0;
}
//----------------------------------------------------------------------------
void vtkITKKnowledgeBasedClassification::SetInput( vtkImageData *input)
{
  if (input)
    {
    this->SetInputConnection(0,input->GetProducerPort());
    }
}

//----------------------------------------------------------------------------
vtkImageData *vtkITKKnowledgeBasedClassification::GetInput()
{
  return vtkImageData::SafeDownCast(this->GetExecutive()->GetInputData(0, 0));
}


//-------------------------------------------------------------------------
int vtkITKKnowledgeBasedClassification:: RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inInfo),
  vtkInformationVector *vtkNotUsed(outInfo))
{
  this->vtkITKPipelineConnection();

  return 0;
}

//-------------------------------------------------------------------------
int vtkITKKnowledgeBasedClassification:: RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inInfo,
  vtkInformationVector *vtkNotUsed(outInfo))
{
  vtkImageData *input = vtkImageData::SafeDownCast(
    inInfo[0]->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));
  this->ImageDataType = input->GetScalarType();
  this->ImageDimension = input->GetDataDimension();

  return 0;
}
//--------------------------------------------------------------------------
// a macro to evaluate an expression for all scalar types
//--------------------------------------------------------------------------

#define vtkTypeCaseMacro(expression) \
      case VTK_FLOAT: { typedef float VTK_TT; expression; } \
        break;

/*#else

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

#endif */

//--------------------------------------------------------------------------
// Template function to Export image from VTK to ITK
//--------------------------------------------------------------------------
template<class T>
static void vtkITK3DConnectVTKITKPipelines(vtkImageExport *vtkImageExporter,
					   itk::Object    *&itkImageImporter,
					   T *dummy )
{
  typedef itk::Image< T , DIMENSIONS> ImageType;

  typename itk::VTKImageImport< ImageType >::Pointer imageSmartPointer
    = itk::VTKImageImport< ImageType >::New();
  imageSmartPointer->Register();
  typename itk::VTKImageImport< ImageType >* imagePointer
    = imageSmartPointer.GetPointer();
  itkImageImporter = imagePointer;

  ::ConnectVTKToITK( vtkImageExporter,
                     imagePointer );
}

//--------------------------------------------------------------------------
// Template function to export Image from ITK to VTK
//--------------------------------------------------------------------------
template<class T>
static void vtkITK3DConnectITKVTKPipelines(vtkImageImport *vtkImageImporter,
					   itk::Object    *&itkImageExporter,
					   T *dummy )
{
  typedef itk::Image< T , DIMENSIONS > VtkOutputImageType;
  typename itk::VTKImageExport< VtkOutputImageType >::Pointer segmentSmartPointer
    = itk::VTKImageExport< VtkOutputImageType >::New();
  segmentSmartPointer->Register();

  typename itk::VTKImageExport< VtkOutputImageType >* segmentPointer
    =  segmentSmartPointer.GetPointer();
  itkImageExporter = segmentPointer;

  ::ConnectITKToVTK( segmentPointer, vtkImageImporter);
}


//-------------------------------------------------------------------
// This Template Function that Performs the Classification Steps in ITK
//-------------------------------------------------------------------
template< class T >
static void vtkITKSegment(
  int numberOfClasses,
  double *userProvidedInitialMeans,
  double *variance,
  itk::Object* itkImageImporter,
  T *dummy)
{
  typedef itk::Image< T , DIMENSIONS> ImageType;
  typedef itk::VTKImageImport< ImageType > VTKImageImportType;
  typedef itk::ScalarImageKmeansImageFilter< ImageType > KMeansFilterType;

  typename ImageType::Pointer image =
    (dynamic_cast<VTKImageImportType*>(itkImageImporter))->GetOutput();
  typename KMeansFilterType::Pointer kmeansFilter = KMeansFilterType::New();

  const unsigned int useNonContiguousLabels = false;
  unsigned int nClasses = numberOfClasses;

  kmeansFilter->SetInput( image );
  kmeansFilter->SetUseNonContiguousLabels( useNonContiguousLabels );

  for( unsigned k=0; k < nClasses; k++ )
    {
    kmeansFilter->AddClassWithInitialMean( userProvidedInitialMeans[k] );
    }

  kmeansFilter->Update();

  typename  KMeansFilterType::ParametersType
    estimatedMeans = kmeansFilter->GetFinalMeans();
  typedef typename KMeansFilterType::OutputImageType KMeansImageType;
  typedef itk::ImageRegionConstIterator< ImageType > ImageDataIteratorType;
  typedef itk::ImageRegionConstIterator< KMeansImageType > KMeansIteratorType;
  typedef itk::Array< double > CovarianceArrayType;

  ImageDataIteratorType itrOnInputImage(image, image->GetRequestedRegion());
  KMeansIteratorType itrOnKMeansImage(
    kmeansFilter->GetOutput(), image->GetRequestedRegion());

  CovarianceArrayType sumsOfSquares( nClasses );
  CovarianceArrayType sums( nClasses );
  CovarianceArrayType estimatedCovariances( nClasses );
  CovarianceArrayType classCount( nClasses );

  // initialize arrays
  for ( unsigned int i = 0; i < nClasses; ++i )
    {
    sumsOfSquares[i] = 0;
    sums[i] = 0;
    classCount[i] = 0;
    }

  // index into image volume using labels. find sumsOfSquares and sums.
  itrOnInputImage.GoToBegin();
  itrOnKMeansImage.GoToBegin();
  int count = 0 ;

  while( !itrOnInputImage.IsAtEnd() )
    {

    sumsOfSquares[(int)itrOnKMeansImage.Get()] =
      sumsOfSquares[(int)itrOnKMeansImage.Get()] +
      itrOnInputImage.Get() * itrOnInputImage.Get();
    sums[(int)itrOnKMeansImage.Get()] =
      sums[(int)itrOnKMeansImage.Get()] +
      itrOnInputImage.Get();
    classCount[(int)itrOnKMeansImage.Get()] =
      classCount[(int)itrOnKMeansImage.Get()] + 1;
    ++itrOnInputImage;
    ++itrOnKMeansImage;
    count++;
    }

 // CALCULATE THE VARIANCE
  for ( unsigned int i = 0; i < nClasses; ++i )
    {
    estimatedCovariances[i] =
      (sumsOfSquares[i] / classCount[i]) -
      ((sums[i] * sums[i]) / (classCount[i] * classCount[i]));
    if ( estimatedCovariances[i] < 0.0000001 )
      {
      estimatedCovariances[i] = 0.0000001;
      };
    }
  for( unsigned k=0; k < nClasses; k++ )
    {
    userProvidedInitialMeans[k] = estimatedMeans[k];
    variance[k] = sqrt(static_cast<double>(estimatedCovariances[k]));
    }
  for(unsigned int i=0; i < nClasses; i++)
    {
    cout<<"Cluster ["<<i<<"] : "<<endl;
    cout<<"       Mean : "<< estimatedMeans[i]<<endl;
    cout<<" Covariance : "<<
          sqrt(static_cast<double>(estimatedCovariances[i]))<<endl;
    }
}

//----------------------------------------------------------------------------
// The main program to set up the segmentation framework and connect
// pieces together and do Classification
//----------------------------------------------------------------------------

int vtkITKKnowledgeBasedClassification::vtkITKPipelineConnection()
{
  void *dummy = 0;

  if (this->GetMTime() > this->ExecuteTime)
    {
    int ret = this->Initialize();
    if (ret == 1)
      {
      return 0;
      }

    this->InvokeEvent(vtkCommand::StartEvent,NULL);
    // reset Abort flag
    this->AbortExecute = 0;
    this->Progress = 0.0;

    ///////////////////////////////////////////////////////
    // Connect vtk and itk pipelines
    int  extent[6];
    vtkImageExport* vtkImageExporter = vtkImageExport::New();
    vtkImageExporter->SetInput( this->GetInput() );
    this->GetInput()->GetExtent(extent);

    switch (this->ImageDataType)
      {
      vtkTypeCaseMacro( vtkITK3DConnectVTKITKPipelines(
                        vtkImageExporter,
                        this->ITKImageImporter,
                        (VTK_TT *)(dummy)) );

      }


    ///////////////////////////////////////////////////////
    // Performs the segmentation

    switch (this->ImageDataType)
      {
      vtkTypeCaseMacro(vtkITKSegment(
			 this->NumberOfClasses,
			 this->Mean,
			 this->Variance,
			 this->ITKImageImporter,
			 (VTK_TT *)dummy) );
      }

    // No need to connect ITK to VTK as this class doesnot have any
    // output
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
// this is not an initialization routine rather a rountine to check
// for correct data type and components
int vtkITKKnowledgeBasedClassification::Initialize()
{
  vtkImageData* image = this->GetInput();
  if ( image == NULL )
    {
    vtkErrorMacro(<<"Image is not present\n");
    return 1;
    }

   return 0;
}
