/*=========================================================================

  Program:   Atamai Classes for VTK
  Module:    vtkITKImageSegmentation3D.h
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
// .NAME vtkITKImageSegmentation3D - VTK class for ITK image segmentation
// .SECTION Description
// This class encapsulates most of ITK's image segmentation
// functionality within a single VTK class.

#ifndef __vtkITKImageSegmentation3D_h
#define __vtkITKImageSegmentation3D_h

#include "vtkITKImageAlgorithm.h"
#include "vtkProcessObject.h"
#include "itkArray.h"
#include "itkObject.h"

#define DIMENSIONS                                        3
#define VTK_SEGMENTATION_ALGORITHM_INVALID               -1
#define INVALID_PARAMETER_VALUE                          -1000.0

//--------------------------------------------------------------------------
// constants for the Segmentation Algorithms

// LEVELSET
#define VTK_LEVEL_SET_CANNY_EDGE                          0
#define VTK_LEVEL_SET_FAST_MARCHING                       1

// WATERSHED TRANSFORM
#define VTK_WATERSHED                                     2

// KNOLEDGEBASED SEGMENTATION FOR GREY WHITE SEGMENTATION
#define VTK_KNOWLEDGEBASED                                3
#define MAX_NUMBER_OF_CLASSES                             8

#define VTK_NUMBER_OF_ALGORITHMS                          4

#define SEGMENTATION_OUTPUT_TYPE                          5

//---------------------------------------------------------------------------
//
class vtkImageData;
class vtkPointSet;
class vtkImageImport;
class vtkAbstractTransform;
class vtkImageStencil;
class vtkImageStencilData;
class vtkMatrix4x4;
class vtkMatrixToHomogeneousTransform;
class vtkImageToImageStencil;

//---------------------------------------------------------------------------
class VTK_EXPORT vtkITKImageSegmentation3D : public vtkITKImageAlgorithm
{
public:

  static vtkITKImageSegmentation3D *New();
  vtkTypeMacro(vtkITKImageSegmentation3D, vtkITKImageAlgorithm);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Returns the number of inputs an algorithm requires if the
  //algorithm has been set, otherwise return -1
  int GetNumberOfInputs();

  // Description:
  // VTK style set/get input methods
  virtual void SetInput(vtkImageData *input);
  virtual vtkImageData* GetInput();

  // Description:
  // Set/Get the image to use as the input image.
  void SetFeatureImage(vtkImageData *input);
  vtkImageData* GetFeatureImage();

  // Description:
  // Set/Get the label image used to initialize the segmentation.
  void SetLabelImage(vtkImageData *input);
  vtkImageData *GetLabelImage();

  // Description:
  // Set/Get whether the algorithm needs a label image.
  int GetNeedsLabelImage(int algorithm);

  // Description:
  // Set/Get a vtkPointSet that cointains the seed points.
  void SetSeedPoints(vtkPointSet *input);
  vtkPointSet *GetSeedPoints();

  // Description:
  // Set/Get whether the algorithm needs seed points.
  int GetNeedsSeedPoints(int algorithm);

  // Description:
  // Get the label image that was produced by the segmentation.
  vtkImageData *GetOutput();

  // Description:
  // Get the stencil that was produced by the segmentation.
  vtkImageStencilData *GetStencilOutput();

  // Description:
  // Set/Get the segmentation Algorithm name to use.
  void SetAlgorithm(int algorithm);
  int GetAlgorithm();

  // Description:
  // Get Number of Algorithms available
  int GetNumberOfAlgorithms() { return 3; };

  // Description:
  // Get the name of the specified Segmentation Algorithm.
  // A null is returned if i+1 is greater than number of available
  // Segmetation  methods.
  const char *GetAlgorithmName(int i);

  // Description:
  // Returns the Number of Parameters for the Algorithm
  int GetNumberOfParameters(int algorithm);

  // Description:
  // Set The default parameter Values
  int SetDefaultParameterValues();

  // Description:
  // Get the parameter name of the specified Algorithm.
  // A null is returned if i+1 is greater than number of available
  // parameters.
  const char* GetParameterName(int algorithm, int parameter);

  // Description:
  // Set the parameter for specific algorithm.
  void SetParameter(const char *, double);
  double GetDefaultParameterValue(int i);

  // Description:
  // Next few functions are for testing the vtkAtamaiAlgorithmGUI
  int GetNumberOfTechniques() { return 4; }
  const char* GetTechniqueName(int techniqueID);

  // Description:
  // Set the threshold to use for stencil generation.
  vtkSetMacro(StencilThreshold, double);
  vtkGetMacro(StencilThreshold, double);

protected:

  vtkITKImageSegmentation3D();
  ~vtkITKImageSegmentation3D();

  // Description:
  // Initialize prior to starting the segmentation.  This will
  // generate an error if any of the components are incompatible.
  // This is automatically called by StartSegmentation.
  int Initialize();

  // Description:
  // this connects the whole pipeline and is called from RequestData()
  virtual int SetVTKITKPipelineConnection();

  // Description:
  // Functions overridden from Superclass
  virtual int RequestData(vtkInformation *request,
                          vtkInformationVector **inInfo,
                          vtkInformationVector *outInfo);
  virtual int RequestInformation(vtkInformation *request,
                                 vtkInformationVector **inInfo,
                                 vtkInformationVector *outInfo);
  virtual int FillInputPortInformation(int port, vtkInformation* info);
  virtual int FillOutputPortInformation(int port, vtkInformation* info);
  virtual void ConvertWorldToVoxel(double w[3], double v[3]);

private:

  vtkITKImageSegmentation3D(const vtkITKImageSegmentation3D&);
  void operator=(const vtkITKImageSegmentation3D&);

  int                               FeatureImageDimension;
  int                               LabelImageDimension;
  int                               ParameterChanged;
  int                               Algorithm;
  int                               Initialized;
  int                               FeatureImageDataType;
  int                               LabelImageDataType;
  double                            StencilThreshold;
  vtkImageImport *                  VTKImageImporter;
  vtkImageToImageStencil *          ImageToStencil;
  vtkTimeStamp                      ExecuteTime;


  //BTX
  itk::Object*                      ITKImageImporter;
  itk::Object*                      ITKLabelImageImporter;
  itk::Object*                      ITKImageExporter;
  itk::Array< double >              Parameters;
  itk::Array< double >              DefaultParameterValues;
  std::vector< itk::Object* >       Filters;
  //ETX
};

#endif //__vtkITKImageSegmentation3D_h
