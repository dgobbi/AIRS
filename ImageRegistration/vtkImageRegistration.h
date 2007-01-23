/*=========================================================================

  Program:   Atamai Classes for VTK
  Module:    $RCSfile: vtkImageRegistration.h,v $

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
// .NAME vtkImageRegistration - Atamai linear registration class
// .SECTION Description
// This class is VTK image registration class.

#ifndef __vtkImageRegistration_h
#define __vtkImageRegistration_h

#include "vtkAlgorithm.h"
#include <vtksys/SystemTools.hxx>

#define VTK_TYPE_INVALID                               -1

// constants for the optimizers
#define VTK_OPTIMIZER_AMOEBA                            0
#define VTK_NUMBER_OF_OPTIMIZERS                        1

// constants for the metrics
#define VTK_METRIC_NORMALIZED_CROSS_CORRELATION         0
#define VTK_METRIC_NORMALIZED_MUTUAL_INFORMATION        1
#define VTK_NUMBER_OF_METRICS                           2

//constants for the interpolators
#define VTK_INTERPOLATOR_NEAREST_NEIGHBOR               0
#define VTK_INTERPOLATOR_LINEAR                         1
#define VTK_INTERPOLATOR_CUBIC                          2
#define VTK_NUMBER_OF_INTERPOLATORS                     3

// constants for the transform types
#define VTK_TRANSFORM_CENTERED                          0
#define VTK_NUMBER_OF_TRANSFORMS                        1

// constants for the preprocessor
#define VTK_PREPROCESSOR_MI                             0
#define VTK_PREPROCESSOR_NC                             1
#define VTK_NUMBER_OF_PREPROCESSORS                     2

class vtkImageData;
class vtkAbstractTransform;
class vtkImageStencilData;
class vtkMatrix4x4;
class vtkMatrixToHomogeneousTransform;
class vtkImageReslice;
class vtkImageGaussianSmooth;
class vtkImageShiftScale;
class vtkTransform;
class vtkImageAccumulate;
class vtkImageRangeCalculator;
class vtkInformation;

class VTK_EXPORT vtkImageRegistration : public vtkAlgorithm
{
public:
  vtkTypeMacro(vtkImageRegistration, vtkAlgorithm);
  static vtkImageRegistration *New();
  void PrintSelf(ostream& os, vtkIndent indent);

  //BTX
#define RegistrationInfoStruct vtkImageRegistration::RegistrationInfo
  class RegistrationInfo
  {
  public:
    vtkImageRegistration            *ImageRegistration;
    vtkObject                       *Optimizer;
    vtkAbstractTransform            *Transform;
    vtkObject                       *Metric;
    vtkImageReslice                 *Interpolator;
    int                              OptimizerType;
    int                              MetricType;
    int                              InterpolatorType;
    int                              TransformType;
    vtkstd::vector<double>*          TransformParametersPointer;
  };
  //ETX

  // Description:
  // Set/Get the image to use as the fixed/target image.
  void SetFixedImage(vtkImageData *input);
  vtkImageData *GetFixedImage();

  // Description:
  // Set/Get the stencil to apply as the FixedImageRegion.
  void SetFixedImageStencil(vtkImageStencilData *stencil);
  vtkImageStencilData *GetFixedImageStencil() ;

  // Description:
  // Set/Get the image to use as the moving/source image.
  void SetMovingImage(vtkImageData *input);
  vtkImageData *GetMovingImage();


  // Description:
  // Set/Get the optimizer to use.
  vtkSetMacro(OptimizerType, int);
  vtkGetMacro(OptimizerType, int);

  // Description:
  // Get the name of the specified optimizer.
  // A null is returned if i+1 is greater than number of available
  // optimization methods.
  const char *GetOptimizerName(int i);

  // Description:
  // Get the parameter name of the specified optimizer.
  // A null is returned if i+1 is greater than number of available
  // parameters.
  const char* GetOptimizerParameterName(int optimizer, int parameter);

  // Description:
  // Set the parameter for specific optimizer.
  void SetOptimizerParameter(const char *, double);

  // Description:
  // Set the image registration metric.
  vtkSetMacro(MetricType, int);
  vtkGetMacro(MetricType, int);

  // Description:
  // Get the name of the specified metric.
  // A null is returned if i+1 is greater than number of available
  // metrics.
  const char *GetMetricName(int i);

  // Description:
  // Get the parameter name of the specified metric.
  // A null is returned if i+1 is greater than number of available
  // parameter.
  const char *GetMetricParameterName(int metric, int parameter);

  // Description:
  // Set the parameter for specific optimizer.
  void SetMetricParameter(const char *, double);

  // Description:
  // Set the image interpolator.
  vtkSetMacro(InterpolatorType, int);
  vtkGetMacro(InterpolatorType, int);

  // Description:
  // Get the name of the specified interpolator.
  // A null is returned if i+1 is greater than number of available
  // interpolators.
  const char *GetInterpolatorName(int i);

  // Description:
  // Set the transform type.
  vtkSetMacro(TransformType, int);
  vtkGetMacro(TransformType, int);
  
  // Description:
  // Get the name of the specified transform type.
  // A null is returned if i+1 is greater than number of available
  // transform types.
  const char *GetTransformName(int i);

  // Description:
  // Get the parameter name of the specified transform type.
  // A null is returned if i+1 is greater than number of available
  // parameters.
  const char *GetTransformParameterName(int transform, int parameter);

  // Description:
  // Set the VTK transform that will be used to initialize the
  // ITK transform.
  // void SetInitialTransform(vtkAbstractTransform *);

  // Description:
  // Set the parameter for specific optimizer.
  void SetTransformParameter(const char *, double);

  // Description:
  // Get the last transform that was produced by the optimizer.
  double GetTransformParameter(const char *);

  // Description:
  // Get the last transform that was produced by the optimizer.
  vtkTransform *GetLastTransform();
  
  // Description:
  // Set the image registration metric.
  vtkSetMacro(PreprocessorType, int);
  vtkGetMacro(PreprocessorType, int);

  // Description:
  // Get the name of the specified metric.
  // A null is returned if i+1 is greater than number of available
  // metrics.
  const char *GetPreprocessorName(int i);

  // Description:
  // Get the parameter name of the specified metric.
  // A null is returned if i+1 is greater than number of available
  // parameter.
  const char *GetPreprocessorParameterName(int metric, int parameter);

  // Description:
  // Set the parameter for specific optimizer.
  void SetPreprocessorParameter(const char *, double);

  // Description:
  // Start registration.  The registration will run to completion,
  // according to the optimization parameters that were set.  To
  // see intermediate results, set an observer for ProgressEvents.
  int UpdateRegistration();

protected:
  vtkImageRegistration();
  ~vtkImageRegistration();

  int  Initialize();
  void InitializeTransform();
  void InitializePreprocessor();
  void InitializeMetric();
  void InitializeOptimizer();
  int  ExecuteRegistration();

  // Functions overridden from Superclass
  virtual int ProcessRequest(vtkInformation *, 
                             vtkInformationVector **,
                             vtkInformationVector *);
  virtual int RequestData(vtkInformation *, 
			  vtkInformationVector **,
			  vtkInformationVector *);
  virtual int RequestUpdateExtent(vtkInformation *vtkNotUsed(request), 
                                 vtkInformationVector **inInfo,
                                 vtkInformationVector *vtkNotUsed(outInfo));
  virtual int RequestInformation(vtkInformation *vtkNotUsed(request), 
                                 vtkInformationVector **inInfo,
                                 vtkInformationVector *vtkNotUsed(outInfo));
  virtual int FillInputPortInformation(int port, vtkInformation* info);
  virtual int FillOutputPortInformation(int port, vtkInformation* info);

  int                              OptimizerType;
  int                              MetricType;
  int                              InterpolatorType;
  int                              TransformType;
  int                              PreprocessorType;

  int                              CurrentIteration;
  double                           Value;
  double                           CurPosition[12];
  double                           SourceImageRescaleIntercept;
  double                           SourceImageRescaleSlope;
  double                           TargetImageRescaleIntercept;
  double                           TargetImageRescaleSlope;

  RegistrationInfo                 RegistrationArgs;
  vtkTimeStamp                     ExecuteTime;

  //BTX
  vtkObject                       *Optimizer;
  vtkAbstractTransform            *Transform;
  vtkObject                       *Metric;
  vtkTransform                    *LastTransform;

  vtkImageReslice                 *SourceReslice;
  vtkImageReslice                 *TargetReslice;
  vtkImageGaussianSmooth          *SourceBlur;
  vtkImageGaussianSmooth          *TargetBlur;
  vtkImageShiftScale              *SourceRescale;
  vtkImageShiftScale              *TargetRescale;
  vtkImageAccumulate              *SourceAccumulate;
  vtkImageAccumulate              *TargetAccumulate;
  vtkImageRangeCalculator         *SourceRange;
  vtkImageRangeCalculator         *TargetRange;
  
  vtkstd::vector< double >         MetricParameters;
  vtkstd::vector< double >         OptimizerParameters;
  vtkstd::vector< double >         TransformParameters;
  vtkstd::vector< double >         PreprocessorParameters;
  //ETX

private:
  // Copy constructor and assigment operator are purposely not implemented
  vtkImageRegistration(const vtkImageRegistration&) {};
  void operator=(const vtkImageRegistration&) {};  
};

#endif //__vtkImageRegistration_h
