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
// .NAME vtkImageRegistration - Perform linear image registration.
// .SECTION Description
// This class will find the transformation that registers the source
// image to the target image.

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
  // The image to use as the source image.  The source voxels define
  // the points at which the target image will be interpolated during
  // the registration.
  void SetSourceImageInputConnection(vtkAlgorithmOutput *input) {
    this->SetInputConnection(0, input); }
  vtkAlgorithmOutput *GetSourceImageInputConnection() {
    return this->GetInputConnection(0, 0); }
  void SetSourceImage(vtkImageData *input);
  vtkImageData *GetSourceImage();

  // Description:
  // The image to use as the target image.  The target image will be
  // resampled at the source voxel locations, after these locations
  // have been passed through the source-to-target transform.  For
  // best results, the target image should be the one with smaller
  // voxel spacing, because it can be interpolated with higher accuracy
  // than an image with larger voxel spacing.
  void SetTargetImageInputConnection(vtkAlgorithmOutput *input) {
    this->SetInputConnection(1, input); }
  vtkAlgorithmOutput *GetTargetImageInputConnection() {
    return this->GetInputConnection(1, 0); }
  void SetTargetImage(vtkImageData *input);
  vtkImageData *GetTargetImage();

  // Description:
  // Set a stencil to apply to the fixed image, to register by using
  // only a portion of the image.  This can only be done for the fixed image.
  void SetSourceImageStencil(vtkImageStencilData *stencil);
  vtkImageStencilData *GetSourceImageStencil();

  // Optimizer types
  enum
  {
    Amoeba
  };

  // Metric types
  enum
  {
    SquaredDifference,
    CrossCorrelation,
    NormalizedCrossCorrelation,
    NeighborhoodCorrelation,
    MutualInformation,
    NormalizedMutualInformation
  };

  // Interpolator types
  enum
  {
    Nearest,
    Linear,
    Cubic,
    Sinc
  };

  // Transform types
  enum
  {
    Translation,
    Rigid,
    Similarity,
    ScaleSourceAxes,
    ScaleTargetAxes,
    Affine
  };

  // Initializer types
  enum
  {
    None,
    Centered
  };

  // Description:
  // Set the image registration metric.  The default is normalized
  // cross correlation.
  vtkSetMacro(MetricType, int);
  void SetMetricTypeToSquaredDifference() {
    this->SetMetricType(SquaredDifference); }
  void SetMetricTypeToCrossCorrelation() {
    this->SetMetricType(CrossCorrelation); }
  void SetMetricTypeToNormalizedCrossCorrelation() {
    this->SetMetricType(NormalizedCrossCorrelation); }
  void SetMetricTypeToNeighborhoodCorrelation() {
    this->SetMetricType(NeighborhoodCorrelation); }
  void SetMetricTypeToMutualInformation() {
    this->SetMetricType(MutualInformation); }
  void SetMetricTypeToNormalizedMutualInformation() {
    this->SetMetricType(NormalizedMutualInformation); }
  vtkGetMacro(MetricType, int);

  // Description:
  // Set the optimizer.  The default is Amoeba (Nelder-Mead Simplex).
  vtkSetMacro(OptimizerType, int);
  void SetOptimizerTypeToAmoeba() {
    this->SetOptimizerType(Amoeba); }
  vtkGetMacro(OptimizerType, int);

  // Description:
  // Set the image interpolator.  The default is Linear.
  vtkSetMacro(InterpolatorType, int);
  void SetInterpolatorTypeToNearest() {
    this->SetInterpolatorType(Nearest); }
  void SetInterpolatorTypeToLinear() {
    this->SetInterpolatorType(Linear); }
  void SetInterpolatorTypeToCubic() {
    this->SetInterpolatorType(Cubic); }
  vtkGetMacro(InterpolatorType, int);

  // Description:
  // Set the transform type.  The default is Rigid.  The Similarity
  // transform type adds a universal scale factor, ScaleSourceAxes
  // allows scaling along all three source image axes, ScaleTargetAxes
  // allows scaling along all three target image axes.
  vtkSetMacro(TransformType, int);
  void SetTransformTypeToRigid() {
    this->SetTransformType(Rigid); }
  void SetTransformTypeToSimilarity() {
    this->SetTransformType(Similarity); }
  void SetTransformTypeToScaleSourceAxes() {
    this->SetTransformType(ScaleSourceAxes); }
  void SetTransformTypeToScaleTargetAxes() {
    this->SetTransformType(ScaleTargetAxes); }
  void SetTransformTypeToAffine() {
    this->SetTransformType(Affine); }
  vtkGetMacro(TransformType, int);

  // Description:
  // Set the transform dimensionality.  The default is 3D.
  // Only 2D and 3D are supported.  A 2D transform will only
  // modify the x and y coordinates.
  vtkSetMacro(TransformDimensionality, int);
  void SetTransformDimensionalityTo2D() {
    this->SetTransformDimensionality(2); }
  void SetTransformDimensionalityTo3D() {
    this->SetTransformDimensionality(3); }
  vtkGetMacro(TransformDimensionality, int);

  // Description:
  // Set the initializer type.  The default is None.  The Centered
  // initializer sets an initial translation that will center the
  // images over each other.
  vtkSetMacro(InitializerType, int);
  void SetInitializerTypeToNone() {
    this->SetInitializerType(None); }
  void SetInitializerTypeToCentered() {
    this->SetInitializerType(Centered); }
  vtkGetMacro(InitializerType, int);

  // Description:
  // Set the size of the joint histogram for mutual information.
  // The default size is 64 by 64.
  vtkSetVector2Macro(JointHistogramSize, int);
  vtkGetVector2Macro(JointHistogramSize, int);

  // Description:
  // Set the ranges of the two axes of the joint histogram.
  // By default, the joint histogram covers the full range of data values
  // present in the image (this default is used whenever the first value
  // in the range is greater than the second value in the range).
  vtkSetVector2Macro(SourceImageRange, double);
  vtkSetVector2Macro(TargetImageRange, double);
  vtkGetVector2Macro(SourceImageRange, double);
  vtkGetVector2Macro(TargetImageRange, double);

  // Description:
  // Initialize the transform.  This will also initialize the
  // NumberOfEvaluations to zero.  If a TransformInitializer is
  // set, then only the rotation part of this matrix will be used,
  // and the initial translation will be set from the initializer.
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
  vtkGetMacro(MetricValue, double);

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

  void ComputeImageRange(vtkImageData *data, vtkImageStencilData *stencil,
                         double range[2]);
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
  int                              InitializerType;
  int                              TransformDimensionality;

  int                              MaximumNumberOfIterations;
  double                           MetricTolerance;
  double                           TransformTolerance;
  double                           MetricValue;

  int                              JointHistogramSize[2];
  double                           SourceImageRange[2];
  double                           TargetImageRange[2];

  vtkTimeStamp                     ExecuteTime;

  vtkObject                       *Optimizer;
  vtkAlgorithm                    *Metric;
  vtkObject                       *Interpolator;
  vtkLinearTransform              *Transform;

  vtkMatrix4x4                    *InitialTransformMatrix;
  vtkImageReslice                 *ImageReslice;
  vtkImageShiftScale              *SourceImageQuantizer;
  vtkImageShiftScale              *TargetImageQuantizer;

  vtkImageRegistrationInfo        *RegistrationInfo;

private:
  // Copy constructor and assigment operator are purposely not implemented
  vtkImageRegistration(const vtkImageRegistration&);
  void operator=(const vtkImageRegistration&);
};

#endif //__vtkImageRegistration_h
