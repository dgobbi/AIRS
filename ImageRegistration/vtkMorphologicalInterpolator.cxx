/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMorphologicalInterpolator.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkMorphologicalInterpolator.h"
#include "vtkImageInterpolatorInternals.h"
#include "vtkImageData.h"
#include "vtkDataArray.h"
#include "vtkObjectFactory.h"

#include "vtkTemplateAliasMacro.h"
// turn off 64-bit ints when templating over all types, because
// they cannot be faithfully represented by doubles
# undef VTK_USE_INT64
# define VTK_USE_INT64 0
# undef VTK_USE_UINT64
# define VTK_USE_UINT64 0

vtkStandardNewMacro(vtkMorphologicalInterpolator);

//----------------------------------------------------------------------------
vtkMorphologicalInterpolator::vtkMorphologicalInterpolator()
{
  this->Operation = VTK_IMIOPERATION_DILATE;
  this->Radius[0] = 0.5;
  this->Radius[1] = 0.5;
  this->Radius[2] = 0.5;
}

//----------------------------------------------------------------------------
vtkMorphologicalInterpolator::~vtkMorphologicalInterpolator()
{
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Operation: "
     << this->GetOperationAsString() << "\n";
  os << indent << "Radius: " << this->Radius[0] << " "
     << this->Radius[1] << " " << this->Radius[2] << "\n";
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::ComputeSupportSize(
  const double vtkNotUsed(matrix)[16], int size[3])
{
  this->InternalUpdate();

  size[0] = 1 + 2*static_cast<int>(this->InternalRadius[0]);
  size[1] = 1 + 2*static_cast<int>(this->InternalRadius[1]);
  size[2] = 1 + 2*static_cast<int>(this->InternalRadius[2]);
}

//----------------------------------------------------------------------------
bool vtkMorphologicalInterpolator::IsSeparable()
{
  return true;
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::SetOperation(int mode)
{
  static int minmode = VTK_IMIOPERATION_DILATE;
  static int maxmode = VTK_IMIOPERATION_ERODE;
  mode = ((mode > minmode) ? mode : minmode);
  mode = ((mode < maxmode) ? mode : maxmode);
  if (this->Operation != mode)
    {
    this->Operation = mode;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
const char *vtkMorphologicalInterpolator::GetOperationAsString()
{
  const char *result = "";

  switch (this->Operation)
    {
    case VTK_IMIOPERATION_DILATE:
      result = "Dilate";
      break;
    case VTK_IMIOPERATION_ERODE:
      result = "Erode";
      break;
    }

  return result;
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::SetRadius(
  double x, double y, double z)
{
  if (this->Radius[0] != x ||
      this->Radius[1] != y ||
      this->Radius[2] != z)
    {
    this->Radius[0] = x;
    this->Radius[1] = y;
    this->Radius[2] = z;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::ComputeInternalRadius(
  const double radius[3])
{
  const double rmin = 1e-17;
  const double rmax = (VTK_IMI_KERNEL_SIZE_MAX-1.0)*0.5;

  // clamp the radius to reasonable values
  for (int i = 0; i < 3; i++)
    {
    double r = radius[i];
    // clamp the radius to a reasonable range
    r = (r >= 0 ? r : 0);
    r = (r <= rmax ? r : rmax);
    this->InternalRadius[i] = r;
    // compute the inverse of the radius
    r = (r >= rmin ? r : rmin); 
    this->InternalRadius[i+3] = 1.0/r;
    }
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::InternalDeepCopy(
  vtkAbstractImageInterpolator *a)
{
  vtkMorphologicalInterpolator *obj =
    vtkMorphologicalInterpolator::SafeDownCast(a);
  if (obj)
    {
    this->SetOperation(obj->Operation);
    this->SetRadius(obj->Radius);
    for (int i = 0; i < 6; i++) // InternalRadius has six elements
      {
      this->InternalRadius[i] = obj->InternalRadius[i];
      }
    }
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::InternalUpdate()
{
  this->ComputeInternalRadius(this->Radius);
  this->InterpolationInfo->InterpolationMode = this->Operation;
  this->InterpolationInfo->ExtraInfo = this->InternalRadius;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//  Interpolation subroutines and associated code
//----------------------------------------------------------------------------

namespace {

//----------------------------------------------------------------------------
template<class F, class T>
struct vtkImageMorphInterpolate
{
  static void General(
    vtkInterpolationInfo *info, const F point[3], F *outPtr);
};

//----------------------------------------------------------------------------
template <class F, class T>
void vtkImageMorphInterpolate<F, T>::General(
  vtkInterpolationInfo *info, const F point[3], F *outPtr)
{
  const T *inPtr = static_cast<const T *>(info->Pointer);
  int *inExt = info->Extent;
  vtkIdType *inInc = info->Increments;
  int numscalars = info->NumberOfComponents;
  int operation = info->InterpolationMode;
  double *radius = static_cast<double *>(info->ExtraInfo);
  double *invRadius = &radius[3];

  // find the closest point and offset from closest point
  double fx, fy, fz;
  int inIdX0 = vtkInterpolationMath::Floor(point[0] + 0.5, fx);
  int inIdY0 = vtkInterpolationMath::Floor(point[1] + 0.5, fy);
  int inIdZ0 = vtkInterpolationMath::Floor(point[2] + 0.5, fz);
  fx -= 0.5;
  fy -= 0.5;
  fz -= 0.5;

  // change arrays into locals
  vtkIdType inIncX = inInc[0];
  vtkIdType inIncY = inInc[1];
  vtkIdType inIncZ = inInc[2];

  int minX = inExt[0];
  int maxX = inExt[1];
  int minY = inExt[2];
  int maxY = inExt[3];
  int minZ = inExt[4];
  int maxZ = inExt[5];

  // arrays for the memory offsets
  vtkIdType factX[VTK_IMI_KERNEL_SIZE_MAX];
  vtkIdType factY[VTK_IMI_KERNEL_SIZE_MAX];
  vtkIdType factZ[VTK_IMI_KERNEL_SIZE_MAX];

  // handle the borders, compute the block size
  int rx = static_cast<int>(radius[0]);
  int ry = static_cast<int>(radius[1]);
  int rz = static_cast<int>(radius[2]);
  int xi = inIdX0 - rx;
  int xm = 2*rx+1;
  int yi = inIdY0 - ry;
  int ym = 2*ry+1;
  int zi = inIdZ0 - rz;
  int zm = 2*rz+1;
  int mm = xm;
  mm = ((mm >= ym) ? mm : ym);
  mm = ((mm >= zm) ? mm : zm);

  switch (info->BorderMode)
    {
    case VTK_IMAGE_BORDER_REPEAT:
      {
      int l = 0;
      int m = mm;
      do
        {
        factX[l] = vtkInterpolationMath::Wrap(xi, minX, maxX)*inIncX;
        factY[l] = vtkInterpolationMath::Wrap(yi, minY, maxY)*inIncY;
        factZ[l] = vtkInterpolationMath::Wrap(zi, minZ, maxZ)*inIncZ;
        l++; xi++; yi++; zi++;
        }
      while (--m);
      }
      break;

    case VTK_IMAGE_BORDER_MIRROR:
      {
      int l = 0;
      int m = mm;
      do
        {
        factX[l] = vtkInterpolationMath::Mirror(xi, minX, maxX)*inIncX;
        factY[l] = vtkInterpolationMath::Mirror(yi, minY, maxY)*inIncY;
        factZ[l] = vtkInterpolationMath::Mirror(zi, minZ, maxZ)*inIncZ;
        l++; xi++; yi++; zi++;
        }
      while (--m);
      }
      break;

    default:
      {
      int l = 0;
      int m = mm;
      do
        {
        factX[l] = vtkInterpolationMath::Clamp(xi, minX, maxX)*inIncX;
        factY[l] = vtkInterpolationMath::Clamp(yi, minY, maxY)*inIncY;
        factZ[l] = vtkInterpolationMath::Clamp(zi, minZ, maxZ)*inIncZ;
        l++; xi++; yi++; zi++;
        }
      while (--m);
      }
      break;
    }

  // the squared distances to each pixel
  F fX[VTK_IMI_KERNEL_SIZE_MAX];
  F fY[VTK_IMI_KERNEL_SIZE_MAX];
  F fZ[VTK_IMI_KERNEL_SIZE_MAX];

  double x = -rx - fx;
  double y = -ry - fy;
  double z = -rz - fz;
  double d;
  int ll = 0;
  do
    {
    d = fabs(x) - 0.5; d = (d < 0 ? 0 : d); d *= invRadius[0]; fX[ll] = d*d;
    d = fabs(y) - 0.5; d = (d < 0 ? 0 : d); d *= invRadius[1]; fY[ll] = d*d;
    d = fabs(z) - 0.5; d = (d < 0 ? 0 : d); d *= invRadius[2]; fZ[ll] = d*d;
    ll++; x++; y++; z++;
    }
  while (--mm);

  // check if only one slice in a particular direction
  int multipleY = (minY != maxY);
  int multipleZ = (minZ != maxZ);

  // the limits to use when doing the interpolation
  int k1 = rz - rz*multipleZ;
  int k2 = rz + rz*multipleZ;
  int j1 = ry - ry*multipleY;
  int j2 = ry + ry*multipleY;

  do // loop over components
    {
    F val = inPtr[factZ[rz] + factY[ry] + factX[rx]];
    int k = k1;
    do // loop over z
      {
      F ifz = fZ[k] - (1.0 + VTK_INTERPOLATE_FLOOR_TOL);
      vtkIdType factz = factZ[k];
      int j = j1;
      do // loop over y
        {
        F ify = fY[j];
        F fzy = ifz + ify;
        vtkIdType factzy = factz + factY[j];
        // loop over x
        const T *tmpPtr = inPtr + factzy;
        const F *tmpfX = fX;
        const vtkIdType *tmpfactX = factX;
        int l = xm;
        if (operation == VTK_IMIOPERATION_DILATE)
          { // dilation uses the max() operator
          do
            {
            if (fzy + tmpfX[0] < 0)
              {
              F tmpval = tmpPtr[tmpfactX[0]];
              val = (val > tmpval ? val : tmpval); // max(val,tmpval)
              }
            tmpfX++;
            tmpfactX++;
            }
          while (--l);
          }
        else // (operation == VTK_IMIOPERATION_ERODE)
          { // erosion uses the min() operator
          do
            {
            if (fzy + tmpfX[0] < 0)
              {
              F tmpval = tmpPtr[tmpfactX[0]];
              val = (val < tmpval ? val : tmpval); // mix(val,tmpval)
              }
            tmpfX++;
            tmpfactX++;
            }
          while (--l);
          }
        }
      while (++j <= j2);
      }
    while (++k <= k2);

    *outPtr++ = val;
    inPtr++;
    }
  while (--numscalars);
}

//----------------------------------------------------------------------------
// Get the interpolation function for the specified data types
template<class F>
void vtkMorphologicalInterpolatorGetInterpolationFunc(
  void (**interpolate)(vtkInterpolationInfo *, const F [3], F *),
  int dataType, int vtkNotUsed(interpolationMode))
{
  switch (dataType)
    {
    vtkTemplateAliasMacro(
      *interpolate =
        &(vtkImageMorphInterpolate<F, VTK_TT>::General)
      );
    default:
      *interpolate = 0;
    }
}

//----------------------------------------------------------------------------
// Interpolation for precomputed weights

template <class F, class T>
struct vtkImageMorphRowInterpolate
{
  static void General(
    vtkInterpolationWeights *weights, int idX, int idY, int idZ,
    F *outPtr, int n);
};

//--------------------------------------------------------------------------
// helper function for high-order interpolation
template<class F, class T>
void vtkImageMorphRowInterpolate<F, T>::General(
  vtkInterpolationWeights *weights, int idX, int idY, int idZ,
  F *outPtr, int n)
{
  int stepX = weights->KernelSize[0];
  int stepY = weights->KernelSize[1];
  int stepZ = weights->KernelSize[2];
  idX *= stepX;
  idY *= stepY;
  idZ *= stepZ;
  int rx = (stepX >> 1);
  int ry = (stepY >> 1);
  int rz = (stepZ >> 1);
  const F *fX = static_cast<F *>(weights->Weights[0]) + idX;
  const F *fY = static_cast<F *>(weights->Weights[1]) + idY;
  const F *fZ = static_cast<F *>(weights->Weights[2]) + idZ;
  const vtkIdType *factX = weights->Positions[0] + idX;
  const vtkIdType *factY = weights->Positions[1] + idY;
  const vtkIdType *factZ = weights->Positions[2] + idZ;
  const T *inPtr = static_cast<const T *>(weights->Pointer);
  int operation = weights->InterpolationMode;

  int numscalars = weights->NumberOfComponents;
  for (int i = n; i > 0; --i)
    {
    const T *inPtr0 = inPtr;
    int c = numscalars;
    do // loop over components
      {
      F val = inPtr[factZ[rz] + factY[ry] + factX[rx]];
      int k = 0;
      do // loop over z
        {
        F ifz = fZ[k] - (1.0 + VTK_INTERPOLATE_FLOOR_TOL);
        vtkIdType factz = factZ[k];
        int j = 0;
        do // loop over y
          {
          F ify = fY[j];
          F fzy = ifz + ify;
          vtkIdType factzy = factz + factY[j];
          // loop over x
          const T *tmpPtr = inPtr0 + factzy;
          const F *tmpfX = fX;
          const vtkIdType *tmpfactX = factX;
          int l = stepX;
          if (operation == VTK_IMIOPERATION_DILATE)
            { // dilation uses the max() operator
            do
              {
              if (fzy + tmpfX[0] < 0)
                {
                F tmpval = tmpPtr[tmpfactX[0]];
                val = (val > tmpval ? val : tmpval); // max(val,tmpval)
                }
              tmpfX++;
              tmpfactX++;
              }
            while (--l);
            }
          else // (operation == VTK_IMIOPERATION_ERODE)
            { // erosion uses the min() operator
            do
              {
              if (fzy + tmpfX[0] < 0)
                {
                F tmpval = tmpPtr[tmpfactX[0]];
                val = (val < tmpval ? val : tmpval); // mix(val,tmpval)
                }
              tmpfX++;
              tmpfactX++;
              }
            while (--l);
            }
          }
        while (++j < stepY);
        }
      while (++k < stepZ);

      *outPtr++ = val;
      inPtr0++;
      }
    while (--c);

    factX += stepX;
    fX += stepX;
    }
}

//----------------------------------------------------------------------------
// get row interpolation function for different interpolation modes
// and different scalar types
template<class F>
void vtkMorphologicalInterpolatorGetRowInterpolationFunc(
  void (**summation)(vtkInterpolationWeights *weights, int idX, int idY,
                     int idZ, F *outPtr, int n),
  int scalarType, int vtkNotUsed(interpolationMode))
{
  switch (scalarType)
    {
    vtkTemplateAliasMacro(
      *summation = &(vtkImageMorphRowInterpolate<F,VTK_TT>::General)
      );
    default:
      *summation = 0;
    }
}

//----------------------------------------------------------------------------
template<class F>
void vtkMorphologicalInterpolatorPrecomputeWeights(
  const F newmat[16], const int outExt[6], int clipExt[6],
  const F bounds[6], vtkInterpolationWeights *weights)
{
  double *radius = static_cast<double *>(weights->ExtraInfo);
  double *invRadius = &radius[3];
  weights->WeightType = vtkTypeTraits<F>::VTKTypeID();

  // set up input positions table for interpolation
  bool validClip = true;
  for (int j = 0; j < 3; j++)
    {
    // set k to the row for which the element in column j is nonzero,
    // and set matrow to the elements of that row
    int k = 0;
    const F *matrow = newmat;
    while (k < 3 && matrow[j] == 0)
      {
      k++;
      matrow += 4;
      }

    // get the extents
    clipExt[2*j] = outExt[2*j];
    clipExt[2*j + 1] = outExt[2*j + 1];
    int minExt = weights->Extent[2*k];
    int maxExt = weights->Extent[2*k + 1];
    F minBounds = bounds[2*k];
    F maxBounds = bounds[2*k + 1];

    int r = static_cast<int>(radius[j]);
    int step = r*2 + 1;
    if (minExt == maxExt)
      {
      step = 1;
      }

    // allocate space for the weights
    vtkIdType size = step*(outExt[2*j+1] - outExt[2*j] + 1);
    vtkIdType *positions = new vtkIdType[size];
    positions -= step*outExt[2*j];
    F *constants = new F[size];
    constants -= step*outExt[2*j];

    weights->KernelSize[j] = step;
    weights->Positions[j] = positions;
    weights->Weights[j] = constants;
    weights->WeightExtent[2*j] = outExt[2*j];
    weights->WeightExtent[2*j+1] = outExt[2*j+1];

    int region = 0;
    for (int i = outExt[2*j]; i <= outExt[2*j+1]; i++)
      {
      F point = matrow[3] + i*matrow[j];

      F f;
      int idx = vtkInterpolationMath::Floor(point + 0.5, f);
      f -= 0.5;
      if (step > 1)
        {
        idx -= r;
        }

      int inId[VTK_IMI_KERNEL_SIZE_MAX];

      int l = 0;
      switch (weights->BorderMode)
        {
        case VTK_IMAGE_BORDER_REPEAT:
          do
            {
            inId[l] = vtkInterpolationMath::Wrap(idx++, minExt, maxExt);
            }
          while (++l < step);
          break;

        case VTK_IMAGE_BORDER_MIRROR:
          do
            {
            inId[l] = vtkInterpolationMath::Mirror(idx++, minExt, maxExt);
            }
          while (++l < step);
          break;

        default:
           do
            {
            inId[l] = vtkInterpolationMath::Clamp(idx++, minExt, maxExt);
            }
          while (++l < step);
          break;
        }

      // compute the weights and offsets
      vtkIdType inInc = weights->Increments[k];
      if (step == 1)
        {
        positions[step*i] = inId[0]*inInc;
        constants[step*i] = 0;
        }
      else
        {
        int ll = 0;
        F x = -f - r;
        do
          {
          F d = fabs(x) - 0.5;
          d = (d < 0 ? 0 : d);
          d *= invRadius[j];
          constants[step*i + ll] = d*d;
          positions[step*i + ll] = inId[ll]*inInc;
          x++;
          }
        while (++ll < step);
        }

      if (point >= minBounds && point <= maxBounds)
        {
        if (region == 0)
          { // entering the input extent
          region = 1;
          clipExt[2*j] = i;
          }
        }
      else
        {
        if (region == 1)
          { // leaving the input extent
          region = 2;
          clipExt[2*j+1] = i - 1;
          }
        }
      }

    if (region == 0 || clipExt[2*j] > clipExt[2*j+1])
      { // never entered input extent!
      validClip = false;
      }
    }

  if (!validClip)
    {
    // output extent doesn't itersect input extent
    for (int j = 0; j < 3; j++)
      {
      clipExt[2*j] = outExt[2*j];
      clipExt[2*j + 1] = outExt[2*j] - 1;
      }
    }
}

//----------------------------------------------------------------------------
} // ends anonymous namespace

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::GetInterpolationFunc(
  void (**func)(vtkInterpolationInfo *, const double [3], double *))
{
  vtkMorphologicalInterpolatorGetInterpolationFunc(
    func, this->InterpolationInfo->ScalarType, this->Operation);
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::GetInterpolationFunc(
  void (**func)(vtkInterpolationInfo *, const float [3], float *))
{
  vtkMorphologicalInterpolatorGetInterpolationFunc(
    func, this->InterpolationInfo->ScalarType, this->Operation);
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::GetRowInterpolationFunc(
  void (**func)(vtkInterpolationWeights *, int, int, int, double *, int))
{
  vtkMorphologicalInterpolatorGetRowInterpolationFunc(
    func, this->InterpolationInfo->ScalarType, this->Operation);
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::GetRowInterpolationFunc(
  void (**func)(vtkInterpolationWeights *, int, int, int, float *, int))
{
  vtkMorphologicalInterpolatorGetRowInterpolationFunc(
    func, this->InterpolationInfo->ScalarType, this->Operation);
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::PrecomputeWeightsForExtent(
  const double matrix[16], const int extent[6], int newExtent[6],
  vtkInterpolationWeights *&weights)
{
  weights = new vtkInterpolationWeights(*this->InterpolationInfo);

  vtkMorphologicalInterpolatorPrecomputeWeights(
    matrix, extent, newExtent, this->StructuredBoundsDouble, weights);
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::PrecomputeWeightsForExtent(
  const float matrix[16], const int extent[6], int newExtent[6],
  vtkInterpolationWeights *&weights)
{
  weights = new vtkInterpolationWeights(*this->InterpolationInfo);

  vtkMorphologicalInterpolatorPrecomputeWeights(
    matrix, extent, newExtent, this->StructuredBoundsFloat, weights);
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::FreePrecomputedWeights(
  vtkInterpolationWeights *&weights)
{
  this->Superclass::FreePrecomputedWeights(weights);
}
