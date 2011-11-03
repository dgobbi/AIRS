/*=========================================================================

  Program:   Atamai Classes for VTK
  Module:    $RCSfile: vtkITKKnowledgeBasedClassification.h,v $
  Creator:   Piali Das <pdas@atamai.com>
  Language:  C++
  Author:    $Author: pdas $
  Date:      $Date: 2006/11/13 18:43:56 $
  Version:   $Revision: 1.1 $

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
// .NAME vtkITKImageClassification3D - VTK class for ITK image segmentation
// .SECTION Description
// This class encapsulates most of ITK's image segmentation
// functionality within a single VTK class.

#ifndef __vtkITKKnowledgeBasedClassification_h
#define __vtkITKKnowledgeBasedClassification_h

#include "vtkITKImageAlgorithm.h"

#include "itkArray.h"
#include "itkObject.h"
#include "vtkProcessObject.h"

class vtkImageData;
class vtkPointSet;
class vtkImageImport;
class vtkAbstractTransform;
class vtkImageStencil;
class vtkMatrix4x4;
class vtkMatrixToHomogeneousTransform;

#define DIMENSIONS 3
//---------------------------------------------------------------------------
//
class vtkITKKnowledgeBasedClassification : public vtkITKImageAlgorithm
{
public:

  static vtkITKKnowledgeBasedClassification *New();
  vtkTypeMacro(vtkITKKnowledgeBasedClassification, vtkProcessObject);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

    // Description:
  // Set/Get the image to use as the input image.
  void SetInput( vtkImageData *input);
  vtkImageData* GetInput( );

  //Desription:
  // Set/Get the mean of the k-th cluster
  void SetMean( int k, double mean );
  double  GetMean( int k );

  void SetVariance( int k, double var );
  double GetVariance(int k);

  void SetNumberOfClasses( const int n);
  // Description:
  // Start Classification.  The segmentation will run to completion,
  // using the selected algorithm and selected corresponding
  // parameters. Intermediate steps cannot be seen.
  int vtkITKPipelineConnection();


 protected:

  vtkITKKnowledgeBasedClassification();
  ~vtkITKKnowledgeBasedClassification();

  // Description:
  // Overridden the superclass RequestData() Method
  virtual int RequestData(vtkInformation *,
			  vtkInformationVector **,
			  vtkInformationVector *);
  virtual int RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inInfo,
  vtkInformationVector *vtkNotUsed(outInfo));
  // Description:
  // Initialize prior to starting the segmentation.  This will
  // generate an error if any of the components are incompatible.
  // This is automatically called by StartClassification.
  int Initialize();

private:

  // Copy constructor and assigment operator are purposely not implemented
  vtkITKKnowledgeBasedClassification(const vtkITKKnowledgeBasedClassification&) {};
  void operator=(const vtkITKKnowledgeBasedClassification&) {};

  int                               ImageDimension;
  int                               NumberOfClasses;
  double                           *Mean;
  double                           *Variance;
  int                               ImageDataType ;
  vtkImageImport *                  VTKImageImporter;
  vtkTimeStamp                      ExecuteTime;

  //BTX
  itk::Object*                      ITKImageImporter;
  itk::Object*                      ITKImageExporter;
  itk::Array< double >              Parameters;
  //ETX
};

#endif //__vtkITKKnowledgeBasedClassification_h
