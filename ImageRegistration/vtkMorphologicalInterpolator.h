/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMorphologicalInterpolator.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkMorphologicalInterpolator - morphological image interpolation
// .SECTION Description
// vtkMorphologicalInterpolator is an interpolator that applies a
// morphological filter (e.g. erosion or dilation) to the image.
// .SECTION Thanks
// Thanks to David Gobbi, Department of Radiology, University of Calgary
// for providing this class.
// .SECTION See also
// vtkImageReslice


#ifndef vtkMorphologicalInterpolator_h
#define vtkMorphologicalInterpolator_h

#include "vtkImageRegistrationModule.h" // For export macro
#include "vtkAbstractImageInterpolator.h"

#define VTK_IMIOPERATION_DILATE 0
#define VTK_IMIOPERATION_ERODE  1
#define VTK_IMI_KERNEL_SIZE_MAX 32

class vtkImageData;
struct vtkInterpolationInfo;

class VTKIMAGEREGISTRATION_EXPORT vtkMorphologicalInterpolator :
  public vtkAbstractImageInterpolator
{
public:
  static vtkMorphologicalInterpolator *New();
  vtkTypeMacro(vtkMorphologicalInterpolator, vtkAbstractImageInterpolator);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the operation to apply within the chosen radius (default: Dilate).
  // If the operation is Dilate, the output pixel is set to the
  // brightest of the input pixels, causing bright regions to be dilated.
  // If the operation is Erode, the output pixel is set to the darkest
  // of the input pixels, causing bright regions to be eroded.
  virtual void SetOperation(int operation);
  void SetOperationToDilate() {
    this->SetOperation(VTK_IMIOPERATION_DILATE); }
  void SetOperationToErode() {
    this->SetOperation(VTK_IMIOPERATION_ERODE); }
  int GetOperation() { return this->Operation; }
  virtual const char *GetOperationAsString();

  // Description:
  // Set the radius of the ellipsoidal structuring element.
  void SetRadius(double x, double y, double z);
  void SetRadius(const double f[3]) {
    this->SetRadius(f[0], f[1], f[2]); }
  void GetRadius(double f[3]) {
    f[0] = this->Radius[0];
    f[1] = this->Radius[1];
    f[2] = this->Radius[2]; }
  double *GetRadius() { return this->Radius; }

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
  vtkMorphologicalInterpolator();
  ~vtkMorphologicalInterpolator();

  // Description:
  // Update the interpolator.
  virtual void InternalUpdate();

  // Description:
  // Copy the interpolator.
  virtual void InternalDeepCopy(vtkAbstractImageInterpolator *obj);

  // Description:
  // Compute the InternalRadius from the user supplied Radius.
  virtual void ComputeInternalRadius(const double radius[3]);

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

  int Operation;
  double Radius[3];
  double InternalRadius[6];

private:
  vtkMorphologicalInterpolator(const vtkMorphologicalInterpolator&);  // Not implemented.
  void operator=(const vtkMorphologicalInterpolator&);  // Not implemented.
};

#endif
