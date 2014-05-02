/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGaussianInterpolator.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkGaussianInterpolator - interpolate with Gaussian kernel
// .SECTION Description
// vtkGaussianInterpolator performs image interpolation with various
// kernels that are derived from the Gaussian function.  The simplest
// kernel is simply the Gaussian itself.  Higher-order kernels use a
// linear combination of the Gaussian and its derivatives to approximate
// a sinc kernel.  These kernels have very favorable filtering
// characteristics (i.e. smooth roll-off in both image space and
// frequency space).  Increasing the order of the kernel allows a better
// approximation of a sinc function, and hence a sharper frequency cutoff.
// These kernels are described in the following publication:
// [1] C. Robert Appledorn, "A New Approach to the Interpolation of Sampled
//     Data," IEEE Transactions on Medical Imaging, 15(3):369-376, 1996.
// .SECTION Thanks
// Thanks to David Gobbi, Department of Radiology, University of Calgary
// for providing this class.
// .SECTION See also
// vtkImageReslice, vtkImageResize


#ifndef __vtkGaussianInterpolator_h
#define __vtkGaussianInterpolator_h

#include "vtkAbstractImageInterpolator.h"

#define VTK_GAUSSIAN_INTERPOLATION 0
#define VTK_APPLEDORN2_INTERPOLATION 1
#define VTK_APPLEDORN6_INTERPOLATION 2
#define VTK_APPLEDORN10_INTERPOLATION 3
#define VTK_GAUSS_KERNEL_SIZE_MAX 32

class vtkImageData;
struct vtkInterpolationInfo;

class VTK_EXPORT vtkGaussianInterpolator :
  public vtkAbstractImageInterpolator
{
public:
  static vtkGaussianInterpolator *New();
  vtkTypeMacro(vtkGaussianInterpolator, vtkAbstractImageInterpolator);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Choose the interpolation kernel.  The default is a simple Gaussian
  // kernel with a standard deviation of 0.399 voxel widths, or exactly
  // 1/sqrt(2*pi).  The other kernels use a linear combination of the
  // derivatives of the Gaussian function to approximate the shape of a
  // sinc kernel.  When the Appledorn10 is used, the RadiusFactors should
  // be set to 4 to properly truncate the kernel.
  virtual void SetKernelType(int itype);
  void SetKernelTypeToGaussian() {
    this->SetKernelType(VTK_GAUSSIAN_INTERPOLATION); }
  void SetKernelTypeToAppledorn2() {
    this->SetKernelType(VTK_APPLEDORN2_INTERPOLATION); }
  void SetKernelTypeToAppledorn6() {
    this->SetKernelType(VTK_APPLEDORN6_INTERPOLATION); }
  void SetKernelTypeToAppledorn10() {
    this->SetKernelType(VTK_APPLEDORN10_INTERPOLATION); }
  int GetKernelType() { return this->KernelType; }
  virtual const char *GetKernelTypeAsString();

  // Description:
  // The radius of the Gaussian kernel will be equal to the product
  // of this RadiusFactor and the BlurFactor.  In other words, increasing
  // the BlurFactor will automatically increase the radius.  The default
  // RadiusFactor is 3, which equates to 7.5 standard deviations.
  // Increasing the RadiusFactor by 1 increases the radius by 2.5
  // standard deviations.  Note that the cut-off function is rectangular,
  // rather than circular.
  void SetRadiusFactors(double x, double y, double z);
  void SetRadiusFactors(const double f[3]) {
    this->SetRadiusFactors(f[0], f[1], f[2]); }
  void GetRadiusFactors(double f[3]) {
    f[0] = this->RadiusFactors[0];
    f[1] = this->RadiusFactors[1];
    f[2] = this->RadiusFactors[2]; }
  double *GetRadiusFactors() { return this->RadiusFactors; }

  // Description:
  // Blur the image by widening the interpolation kernel by the provided
  // BlurFactor.  This can be used to set the standard deviation of the
  // Gaussian, which will be equal to BlurFactor/sqrt(2*pi).  When the
  // BlurFactor is 1, the standard deviation is 0.399, which provides
  // a Gaussian of unit height and unit area.  For higher-order kernels,
  // the frequency cutoff is reduced by the provided BlurFactor.
  // If you turn Antialiasing on, then the blur factors will be computed
  // automatically so that the frequency cutoff is the Nyquist frequency
  // for the output sampling rate.  Blurring increases the computation
  // time because the kernel size increases by the blur factor.
  void SetBlurFactors(double x, double y, double z);
  void SetBlurFactors(const double f[3]) {
    this->SetBlurFactors(f[0], f[1], f[2]); }
  void GetBlurFactors(double f[3]) {
    f[0] = this->BlurFactors[0];
    f[1] = this->BlurFactors[1];
    f[2] = this->BlurFactors[2]; }
  double *GetBlurFactors() { return this->BlurFactors; }

  // Description:
  // Turn on antialiasing.  If antialiasing is on, then the BlurFactors
  // will be computed automatically from the output sampling rate such that
  // that the image will be bandlimited to the Nyquist frequency.  This
  // is only applicable when the interpolator is being used by a resampling
  // filter like vtkImageReslice.  Such a filter will indicate the output
  // sampling by calling the interpolator's ComputeSupportSize() method,
  // which will compute the blur factors at the same time that it computes
  // the support size.
  void SetAntialiasing(int antialiasing);
  void AntialiasingOn() { this->SetAntialiasing(1); }
  void AntialiasingOff() { this->SetAntialiasing(0); }
  int GetAntialiasing() { return this->Antialiasing; }

  // Description:
  // Get the support size for use in computing update extents.  If the data
  // will be sampled on a regular grid, then pass a matrix describing the
  // structured coordinate transformation between the output and the input.
  // Otherwise, pass NULL as the matrix to retrieve the full kernel size.
  virtual void ComputeSupportSize(const double matrix[16], int support[3]);

  // Description:
  // Turn on renormalization (default: off).  Renormalization ensures that
  // the weights computed from the kernel always sum to one, which is usually
  // desirable.  However, for the Gaussian kernel, renormalization creates
  // large negative side-lobes in the frequency response.
  void SetRenormalization(int antialiasing);
  void RenormalizationOn() { this->SetRenormalization(1); }
  void RenormalizationOff() { this->SetRenormalization(0); }
  int GetRenormalization() { return this->Renormalization; }

  // Description:
  // Turn on label interpolation (default: off).  This is used for
  // interpolation of label images.  The interpolation weights for
  // each label are computed, and the label with the greatest weight
  // is chosen as the label for the output voxel.
  void SetLabelling(int lablels);
  void LabellingOn() { this->SetLabelling(1); }
  void LabellingOff() { this->SetLabelling(0); }
  int GetLabelling() { return this->Labelling; }

  // Description:
  // Returns true if the interpolator supports weight precomputation.
  // This will always return true for this interpolator.
  virtual bool IsSeparable();

  // Description:
  // If the data is going to be sampled on a regular grid, then the
  // interpolation weights can be precomputed.  A matrix must be
  // supplied that provides a transformation between the provided
  // extent and the structured coordinates of the input.  This
  // matrix must perform only permutations, scales, and translation,
  // i.e. each of the three columns must have only one non-zero value.
  // A new extent is provided for out-of-bounds checks.
  // THIS METHOD IS THREAD SAFE.
  virtual void PrecomputeWeightsForExtent(
    const double matrix[16], const int extent[6], int newExtent[6],
    vtkInterpolationWeights *&weights);
  virtual void PrecomputeWeightsForExtent(
    const float matrix[16], const int extent[6], int newExtent[6],
    vtkInterpolationWeights *&weights);

  // Description:
  // Free the precomputed weights.  THIS METHOD IS THREAD SAFE.
  virtual void FreePrecomputedWeights(vtkInterpolationWeights *&weights);

protected:
  vtkGaussianInterpolator();
  ~vtkGaussianInterpolator();

  // Description:
  // Update the interpolator.
  virtual void InternalUpdate();

  // Description:
  // Copy the interpolator.
  virtual void InternalDeepCopy(vtkAbstractImageInterpolator *obj);

  // Description:
  // Get the interpolation functions.
  virtual void GetInterpolationFunc(
    void (**doublefunc)(
      vtkInterpolationInfo *, const double [3], double *));
  virtual void GetInterpolationFunc(
    void (**floatfunc)(
      vtkInterpolationInfo *, const float [3], float *));

  // Description:
  // Get the row interpolation functions.
  virtual void GetRowInterpolationFunc(
    void (**doublefunc)(
      vtkInterpolationWeights *, int, int, int, double *, int));
  virtual void GetRowInterpolationFunc(
    void (**floatfunc)(
      vtkInterpolationWeights *, int, int, int, float *, int));

  // Description:
  // Build the lookup tables used for the interpolation.
  virtual void BuildKernelLookupTable();

  // Description:
  // Free the kernel lookup tables.
  virtual void FreeKernelLookupTable();

  int KernelType;
  float *KernelLookupTable[3];
  int KernelSize[3];
  int Antialiasing;
  int Renormalization;
  int Labelling;
  double RadiusFactors[3];
  double BlurFactors[3];
  double LastBlurFactors[3];

private:
  vtkGaussianInterpolator(const vtkGaussianInterpolator&);  // Not implemented.
  void operator=(const vtkGaussianInterpolator&);  // Not implemented.
};

#endif
