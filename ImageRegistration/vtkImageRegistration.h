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

class vtkImageData;
class vtkImageStencilData;
class vtkLinearTransform;
class vtkMatrix4x4;
class vtkImageReslice;
class vtkImageShiftScale;

struct vtkImageRegistrationInfo;

class VTK_EXPORT vtkImageRegistration : public vtkAlgorithm
{
public:
  vtkTypeMacro(vtkImageRegistration, vtkAlgorithm);
  static vtkImageRegistration *New();
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // The image to use as the moving/source image.
  void SetSourceImageInputConnection(vtkAlgorithmOutput *input) {
    this->SetInputConnection(0, input); }
  vtkAlgorithmOutput *GetSourceImageInputConnection() {
    return this->GetInputConnection(0, 0); }
  void SetSourceImage(vtkImageData *input);
  vtkImageData *GetSourceImage();

  // Description:
  // The image to use as the fixed/target image.
  void SetTargetImageInputConnection(vtkAlgorithmOutput *input) {
    this->SetInputConnection(1, input); }
  vtkAlgorithmOutput *GetTargetImageInputConnection() {
    return this->GetInputConnection(1, 0); }
  void SetTargetImage(vtkImageData *input);
  vtkImageData *GetTargetImage();

  // Description:
  // Set a stencil to apply to the fixed image, to register by using
  // only a portion of the image.  This can only be done for the fixed image.
  void SetTargetImageStencil(vtkImageStencilData *stencil);
  vtkImageStencilData *GetTargetImageStencil();

  // Optimizer types
  enum
  {
    Amoeba,
  };

  // Metric types
  enum
  {
    CrossCorrelation,
    NormalizedCrossCorrelation,
    MutualInformation,
    NormalizedMutualInformation,
  };

  // Interpolator types
  enum
  {
    Nearest,
    Linear,
    Cubic,
  };

  // Transform types
  enum
  {
    Rigid,
    Similarity,
  };

  // Description:
  // Set the image registration metric.
  vtkSetMacro(MetricType, int);
  void SetMetricTypeToCrossCorrelation() {
    this->SetMetricType(CrossCorrelation); }
  void SetMetricTypeToNormalizedCrossCorrelation() {
    this->SetMetricType(NormalizedCrossCorrelation); }
  void SetMetricTypeToMutualInformation() {
    this->SetMetricType(MutualInformation); }
  void SetMetricTypeToNormalizedMutualInformation() {
    this->SetMetricType(NormalizedMutualInformation); }
  vtkGetMacro(MetricType, int);

  // Description:
  // Set the optimizer.
  vtkSetMacro(OptimizerType, int);
  void SetOptimizerTypeToAmoeba() {
    this->SetOptimizerType(Amoeba); }
  vtkGetMacro(OptimizerType, int);

  // Description:
  // Set the image interpolator.
  vtkSetMacro(InterpolatorType, int);
  void SetInterpolatorTypeToNearest() {
    this->SetInterpolatorType(Nearest); }
  void SetInterpolatorTypeToLinear() {
    this->SetInterpolatorType(Linear); }
  void SetInterpolatorTypeToCubic() {
    this->SetInterpolatorType(Cubic); }
  vtkGetMacro(InterpolatorType, int);

  // Description:
  // Set the transform type.
  vtkSetMacro(TransformType, int);
  void SetTransformTypeToRigid() {
    this->SetTransformType(Rigid); }
  void SetTransformTypeToSimilarity() {
    this->SetTransformType(Similarity); }
  vtkGetMacro(TransformType, int);

  // Description:
  // Set the size of the joint histogram for mutual information.
  vtkSetVector2Macro(JointHistogramSize, int);
  vtkGetVector2Macro(JointHistogramSize, int);

  // Description:
  // Initialize the transform.  This will also initialize the
  // NumberOfEvaluations to zero.
  void Initialize(vtkMatrix4x4 *matrix);

  // Description:
  // Set the tolerance that the optimizer will apply to the value returned
  // by the image similarity metric.  The default value is 1e-4.
  vtkSetMacro(MetricTolerance, double);
  vtkGetMacro(MetricTolerance, double);

  // Description:
  // Set the tolerance that the optimizer will apply to the transform
  // parameters.  Provide the value to use for the translation parameters,
  // it will automatically be scaled when applied to the rotation and
  // scale parameters.  The default value is 0.1.
  vtkSetMacro(TransformTolerance, double);
  vtkGetMacro(TransformTolerance, double);

  // Description:
  // Set the maximum number of iterations to perform.  The number of metric
  // evaluations per iteration will depend on the optimizer.
  vtkSetMacro(MaximumNumberOfIterations, int);
  vtkGetMacro(MaximumNumberOfIterations, int);

  // Description:
  // Get the number of times that the metric has been evaluated.
  int GetNumberOfEvaluations();

  // Description:
  // Get the last transform that was produced by the optimizer.
  vtkLinearTransform *GetTransform() { return this->Transform; }

  // Description:
  // Get the value that is being minimized.
  vtkGetMacro(Value, double);

  // Description:
  // Iterate the registration.  Returns zero if the termination condition has
  // been reached.
  int Iterate();

  // Description:
  // Start registration.  The registration will run to completion,
  // according to the optimization parameters that were set.  To
  // see intermediate results, set an observer for ProgressEvents.
  // There will be one ProgressEvent per iteration.  This will return
  // zero if the maximum number of iterations was reached before convergence.
  int UpdateRegistration();

protected:
  vtkImageRegistration();
  ~vtkImageRegistration();

  void ComputeImageRange(vtkImageData *data, double range[2]);
  int ExecuteRegistration();

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

  int                              JointHistogramSize[2];
  int                              MaximumNumberOfIterations;
  double                           MetricTolerance;
  double                           TransformTolerance;
  double                           Value;

  vtkTimeStamp                     ExecuteTime;

  vtkObject                       *Optimizer;
  vtkAlgorithm                    *Metric;
  vtkObject                       *Interpolator;
  vtkLinearTransform              *Transform;

  vtkImageReslice                 *ImageReslice;
  vtkImageShiftScale              *SourceImageQuantizer;
  vtkImageShiftScale              *TargetImageQuantizer;

  vtkImageRegistrationInfo        *RegistrationInfo;

private:
  // Copy constructor and assigment operator are purposely not implemented
  vtkImageRegistration(const vtkImageRegistration&) {};
  void operator=(const vtkImageRegistration&) {};
};

#endif //__vtkImageRegistration_h
