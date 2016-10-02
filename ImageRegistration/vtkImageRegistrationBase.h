/*=========================================================================

  Module: vtkImageRegistrationBase.h

  Copyright (c) 2006 Atamai, Inc.
  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageRegistrationBase - Base class for registration classes
// .SECTION Description
// This class provides the basic interface for image registration.

#ifndef vtkImageRegistrationBase_h
#define vtkImageRegistrationBase_h

#include "vtkAlgorithm.h"

class vtkImageData;
class vtkImageStencilData;
class vtkLinearTransform;
class vtkTransform;
class vtkMatrix4x4;
class vtkDoubleArray;
class vtkAbstractImageInterpolator;
class vtkFunctionMinimizer;
class vtkImageSimilarityMetric;

class VTK_EXPORT vtkImageRegistrationBase : public vtkAlgorithm
{
public:
  vtkTypeMacro(vtkImageRegistrationBase, vtkAlgorithm);
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
  // Set a stencil to apply to the fixed image.
  // The images will only be compared within this region.
  void SetSourceImageStencil(vtkImageStencilData *stencil);
  vtkImageStencilData *GetSourceImageStencil();

  // Description:
  // Set a stencil to apply to the moving image.
  // This mask moves with the image.  It is resampled and then combined
  // with the fixed mask in order to generate the final mask.
  void SetTargetImageStencil(vtkImageStencilData *stencil);
  vtkImageStencilData *GetTargetImageStencil();

  // Optimizer types
  enum
  {
    Amoeba,
    Powell
  };

  // Metric types
  enum
  {
    SquaredDifference,
    CrossCorrelation,
    NormalizedCrossCorrelation,
    NeighborhoodCorrelation,
    CorrelationRatio,
    MutualInformation,
    NormalizedMutualInformation
  };

  // Interpolator types
  enum
  {
    Nearest,
    Linear,
    Cubic,
    BSpline,
    Sinc,
    ASinc,
    Label
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
  // Set the image registration metric.  The default is mutual information.
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
  // Set the optimizer.  The default is Powell.
  vtkSetMacro(OptimizerType, int);
  void SetOptimizerTypeToAmoeba() {
    this->SetOptimizerType(Amoeba); }
  void SetOptimizerTypeToPowell() {
    this->SetOptimizerType(Powell); }
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
  void SetInterpolatorTypeToBSpline() {
    this->SetInterpolatorType(BSpline); }
  void SetInterpolatorTypeToSinc() {
    this->SetInterpolatorType(Sinc); }
  void SetInterpolatorTypeToASinc() {
    this->SetInterpolatorType(ASinc); }
  void SetInterpolatorTypeToLabel() {
    this->SetInterpolatorType(Label); }
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
  // Set the tolerance that the optimizer will apply to the cost value.
  // This is a relative, rather than absolute, tolerance.  The default
  // value is 1e-4.
  vtkSetMacro(CostTolerance, double);
  vtkGetMacro(CostTolerance, double);

  // Description:
  // Set the tolerance that the optimizer will apply to the transform
  // parameters.  Provide the value to use for the translation parameters,
  // it will automatically be scaled when applied to the rotation and
  // scale parameters.  The default value is 0.1.
  vtkSetMacro(TransformTolerance, double);
  vtkGetMacro(TransformTolerance, double);

  // Description:
  // Set the maximum number of iterations to perform. Default: 500.
  // The number of metric evaluations per iteration will depend on the
  // optimizer that is used.
  vtkSetMacro(MaximumNumberOfIterations, int);
  vtkGetMacro(MaximumNumberOfIterations, int);

  // Description:
  // Set the maximum number of metric evaluations. Default: 5000.
  // This is usually a more useful limit than the number of iterations.
  // Note the registration will continue until the end of the current
  // iteration.
  vtkSetMacro(MaximumNumberOfEvaluations, int);
  vtkGetMacro(MaximumNumberOfEvaluations, int);

  // Description:
  // Copy all settings from another image registration object.
  void CopySettings(vtkImageRegistrationBase *other);

  // Description:
  // Get the number of times that the metric has been evaluated.
  int GetNumberOfEvaluations();

  // Description:
  // Get the last transform that was produced by the optimizer.
  vtkLinearTransform *GetTransform();

  // Description:
  // Get the current cost value for the registration.
  vtkGetMacro(CostValue, double);

  // Description:
  // Turn this on to collect diagnositic values during registration.
  // This will cause the MetricValues, CostValues, and ParameterValues
  // to be collected during the registration.
  vtkGetMacro(CollectValues, bool);
  vtkSetMacro(CollectValues, bool);
  vtkBooleanMacro(CollectValues, bool);

  // Description:
  // Get an array of all metric values since registration started.
  // There will be one metric value stored for each function evaluation
  // that was performed.
  vtkGetObjectMacro(MetricValues, vtkDoubleArray)

  // Description:
  // Get an array of all cost values since registration started.
  // There will be one value stored for each function evaluation
  // that was performed.
  vtkGetObjectMacro(CostValues, vtkDoubleArray)

  // Description:
  // Get an array of all parameter values since registration started.
  // There will be one set of parameters stored for each function evaluation
  // that was performed.  The function parameters are, in order, the
  // translation parameters, the rotation parameters, the scale followed
  // by the scale anisotropy, and finally the rotation parameters for the
  // axes of the scale parameters.
  vtkGetObjectMacro(ParameterValues, vtkDoubleArray)

  // Description:
  // Initialize the transform.  This will also initialize the
  // NumberOfEvaluations to zero.  If a TransformInitializer is
  // set, then only the rotation part of this matrix will be used,
  // and the initial translation will be set from the initializer.
  virtual void Initialize(vtkMatrix4x4 *matrix) = 0;

  // Description:
  // Iterate the registration.  Returns zero if the termination condition has
  // been reached.
  virtual int Iterate() = 0;

  // Description:
  // Start registration.  The registration will run to completion,
  // according to the optimization parameters that were set.  To
  // see intermediate results, set an observer for ProgressEvents.
  // There will be one ProgressEvent per iteration.  This will return
  // zero if the maximum number of iterations was reached before convergence.
  int UpdateRegistration();

protected:
  vtkImageRegistrationBase();
  ~vtkImageRegistrationBase();

  void ComputeImageRange(vtkImageData *data, vtkImageStencilData *stencil,
                         double range[2]);

  virtual int ExecuteRegistration() = 0;

  // Functions overridden from Superclass
  virtual int ProcessRequest(vtkInformation *,
                             vtkInformationVector **,
                             vtkInformationVector *);
  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *);
  virtual int RequestUpdateExtent(vtkInformation *,
                                 vtkInformationVector **,
                                 vtkInformationVector *);
  virtual int RequestInformation(vtkInformation *,
                                 vtkInformationVector **,
                                 vtkInformationVector *);
  virtual int FillInputPortInformation(int port, vtkInformation* info);
  virtual int FillOutputPortInformation(int port, vtkInformation* info);

  int                              OptimizerType;
  int                              MetricType;
  int                              InterpolatorType;
  int                              TransformType;
  int                              InitializerType;
  int                              TransformDimensionality;

  int                              MaximumNumberOfIterations;
  int                              MaximumNumberOfEvaluations;
  int                              NumberOfEvaluations;
  double                           CostTolerance;
  double                           TransformTolerance;
  double                           CostValue;

  int                              JointHistogramSize[2];
  double                           SourceImageRange[2];
  double                           TargetImageRange[2];

  vtkTimeStamp                     ExecuteTime;

  vtkFunctionMinimizer            *Optimizer;
  vtkImageSimilarityMetric        *Metric;
  vtkAbstractImageInterpolator    *Interpolator;
  vtkTransform                    *Transform;

  bool                             CollectValues;
  vtkDoubleArray                  *MetricValues;
  vtkDoubleArray                  *CostValues;
  vtkDoubleArray                  *ParameterValues;

private:
  // Copy constructor and assigment operator are purposely not implemented
  vtkImageRegistrationBase(const vtkImageRegistrationBase&);
  void operator=(const vtkImageRegistrationBase&);
};

#endif /* vtkImageRegistrationBase_h */
