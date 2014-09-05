/*=========================================================================
  Program:   Atamai Image Registration and Segmentation
  Module:    vtkLabelInterpolator.h

  Copyright (c) 2014 David Gobbi
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

  * Neither the name of the Calgary Image Processing and Analysis Centre
    (CIPAC), the University of Calgary, nor the names of any authors nor
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=========================================================================*/
// .NAME vtkLabelInterpolator - interpolator for label images
// .SECTION Description
// vtkLabelInterpolator is an image interpolator for label images, that is,
// it is meant to be used on images that consist of discrete label values
// rather than continuous intensity values.  It works by computing, for each
// voxel in the output image, the label value that is most probable at that
// voxel location.  It assumes that each input voxel has a Gaussian
// probability distribution.  That is, for each output voxel, the probability
// of the voxel having a specific label value is the sum of the Gaussian
// distributions of all nearby input voxels with that label value.
// .SECTION Thanks
// This class was written by David Gobbi, Calgary Image Procesing and
// Analysis Centre, University of Calgary.
// .SECTION See also
// vtkImageReslice, vtkImageResize


#ifndef __vtkLabelInterpolator_h
#define __vtkLabelInterpolator_h

#include "vtkAbstractImageInterpolator.h"

#define VTK_LABEL_KERNEL_SIZE_MAX 32
#define VTK_LABEL_KERNEL_LABELS_MAX 32

class vtkImageData;
struct vtkInterpolationInfo;

class VTK_EXPORT vtkLabelInterpolator :
  public vtkAbstractImageInterpolator
{
public:
  static vtkLabelInterpolator *New();
  vtkTypeMacro(vtkLabelInterpolator, vtkAbstractImageInterpolator);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

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
  // a Gaussian of unit height and unit area.  If you turn Antialiasing on,
  // then the blur factors will be computed automatically so that the
  // frequency cutoff is the Nyquist frequency for the output sampling rate.
  // Blurring increases the computation time because the kernel size
  // increases by the blur factor.
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
  vtkLabelInterpolator();
  ~vtkLabelInterpolator();

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

  float *KernelLookupTable[3];
  int KernelSize[3];
  int Antialiasing;
  double RadiusFactors[3];
  double BlurFactors[3];
  double LastBlurFactors[3];

private:
  vtkLabelInterpolator(const vtkLabelInterpolator&);  // Not implemented.
  void operator=(const vtkLabelInterpolator&);  // Not implemented.
};

#endif
