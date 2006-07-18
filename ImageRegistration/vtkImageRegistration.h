/*=========================================================================

  Program:   Atamai Classes for VTK
  Module:    $RCSfile: vtkImageRegistration.h,v $
  Creator:   Kevin Wang <kwang@atamai.com>
  Language:  C++
  Author:    $Author: kwang $
  Date:      $Date: 2006/07/18 15:25:40 $
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
// .NAME vtkImageRegistration - VTK class for image registration
// .SECTION Description
// This class is VTK image registration class.

#ifndef __vtkImageRegistration_h
#define __vtkImageRegistration_h

#include "vtkProcessObject.h"
#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION <= 4)
#include <vector>
#else
#include <vtksys/SystemTools.hxx>
#endif

#define VTK_TYPE_INVALID                               -1

//--------------------------------------------------------------------------
// constants for the optimizers
#define VTK_OPTIMIZER_AMOEBA                            0
#define VTK_NUMBER_OF_OPTIMIZERS                        1

//---------------------------------------------------------------------------
// constants for the metrics
#define VTK_METRIC_NORMALIZED_CROSS_CORRELATION         0
#define VTK_METRIC_NORMALIZED_MUTUAL_INFORMATION        1
#define VTK_NUMBER_OF_METRICS                           2

//--------------------------------------------------------------------------
//constants for the interpolators
#define VTK_INTERPOLATOR_NEAREST_NEIGHBOR               0
#define VTK_INTERPOLATOR_LINEAR                         1
#define VTK_INTERPOLATOR_CUBIC                          2
#define VTK_NUMBER_OF_INTERPOLATORS                     3

//---------------------------------------------------------------------------
// constants for the transform types
#define VTK_TRANSFORM_CENTERED                          0
#define VTK_NUMBER_OF_TRANSFORMS                        1

//---------------------------------------------------------------------------
//BTX
class vtkImageData;
class vtkAbstractTransform;
class vtkImageStencilData;
class vtkMatrix4x4;
class vtkMatrixToHomogeneousTransform;
class vtkImageReslice;
//ETX

//---------------------------------------------------------------------------
//
class VTK_EXPORT vtkImageRegistration : public vtkProcessObject
{
public:
  vtkTypeMacro(vtkImageRegistration, vtkProcessObject);
  static vtkImageRegistration *New();
  void PrintSelf(ostream& os, vtkIndent indent);

  //BTX
#define RegistrationInfoStruct vtkImageRegistration::RegistrationInfo
  class RegistrationInfo
  {
  public:
    vtkImageRegistration*            ImageRegistration;
    vtkObject*                       Optimizer;
    vtkAbstractTransform*            Transform;
    vtkObject*                       Metric;
    vtkImageReslice*                 Interpolator;
    int                              OptimizerType;
    int                              MetricType;
    int                              InterpolatorType;
    int                              TransformType;
#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION < 4)
    std::vector<double>*          TransformParametersPointer;
#else
    vtkstd::vector<double>*          TransformParametersPointer;
#endif
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
  vtkAbstractTransform *GetLastTransform();
  
  // Description:
  // Start registration.  The registration will run to completion,
  // according to the optimization parameters that were set.  To
  // see intermediate results, set an observer for ProgressEvents.
  int UpdateRegistration();

protected:

  vtkImageRegistration();
  ~vtkImageRegistration();

  void InitializeTransformAndParameters();
  void InitializeInterpolator();
  void InitializeMetricAndParameters();
  void InitializeOptimizerAndParameters();
  int  ExecuteRegistration();
  int  Initialize();


private:

  // Copy constructor and assigment operator are purposely not implemented
  vtkImageRegistration(const vtkImageRegistration&) {};
  void operator=(const vtkImageRegistration&) {};  

  // 
  int                              OptimizerType;
  int                              MetricType;
  int                              InterpolatorType;
  int                              TransformType;

  RegistrationInfo                 RegistrationArgs;
  
  vtkTimeStamp                     ExecuteTime;

  int                              CurrentIteration;
  double                           Value;
  double                           CurPosition[12];

  //BTX
  vtkObject*                       Optimizer;
  vtkAbstractTransform*            Transform;
  vtkObject*                       Metric;
  vtkImageReslice*                 Interpolator;
  vtkMatrixToHomogeneousTransform* LastTransform;

#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION < 4)
  std::vector< double >            MetricParameters;
  std::vector< double >            OptimizerParameters;
  std::vector< double >            TransformParameters;
#else
  vtkstd::vector< double >         MetricParameters;
  vtkstd::vector< double >         OptimizerParameters;
  vtkstd::vector< double >         TransformParameters;
#endif
  //ETX

};

#endif //__vtkImageRegistration_h
