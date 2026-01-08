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

#include "vtkLabelInterpolator.h"
#include "vtkImageData.h"
#include "vtkDataArray.h"
#include "vtkObjectFactory.h"
#include "vtkVersion.h"

// VTK 9.6 no longer requires use of "Internal" interpolation header
#if VTK_MAJOR_VERSION < 9 || (VTK_MAJOR_VERSION == 9 && VTK_MINOR_VERSION < 6)
#include "vtkImageInterpolatorInternals.h"
#else
#include "vtkInterpolationMath.h"
#endif

#include "vtkTemplateAliasMacro.h"
// turn off 64-bit ints when templating over all types, because
// they cannot be faithfully represented by doubles
# undef VTK_USE_INT64
# define VTK_USE_INT64 0
# undef VTK_USE_UINT64
# define VTK_USE_UINT64 0

// masks for storing window and size in a single integer
#define VTK_INTERPOLATION_WINDOW_MASK        0x0000007f
#define VTK_INTERPOLATION_WINDOW_XBLUR_MASK  0x00008000
#define VTK_INTERPOLATION_WINDOW_XSIZE_MASK  0x00007f00
#define VTK_INTERPOLATION_WINDOW_XSIZE_SHIFT 8
#define VTK_INTERPOLATION_WINDOW_YBLUR_MASK  0x00800000
#define VTK_INTERPOLATION_WINDOW_YSIZE_MASK  0x007f0000
#define VTK_INTERPOLATION_WINDOW_YSIZE_SHIFT 16
#define VTK_INTERPOLATION_WINDOW_ZBLUR_MASK  0x80000000
#define VTK_INTERPOLATION_WINDOW_ZSIZE_MASK  0x7f000000
#define VTK_INTERPOLATION_WINDOW_ZSIZE_SHIFT 24

// kernel lookup table size must be 256*n where n is kernel half-width
// in order to provide sufficient precision for 16-bit images
#define VTK_LABEL_KERNEL_TABLE_DIVISIONS 256

vtkStandardNewMacro(vtkLabelInterpolator);

//----------------------------------------------------------------------------
vtkLabelInterpolator::vtkLabelInterpolator()
{
  this->RadiusFactors[0] = 3;
  this->RadiusFactors[1] = 3;
  this->RadiusFactors[2] = 3;
  this->KernelLookupTable[0] = NULL;
  this->KernelLookupTable[1] = NULL;
  this->KernelLookupTable[2] = NULL;
  this->KernelSize[0] = 6;
  this->KernelSize[1] = 6;
  this->KernelSize[2] = 6;
  this->Antialiasing = 0;
  this->BlurFactors[0] = 1.0;
  this->BlurFactors[1] = 1.0;
  this->BlurFactors[2] = 1.0;
  this->LastBlurFactors[0] = 1.0;
  this->LastBlurFactors[1] = 1.0;
  this->LastBlurFactors[2] = 1.0;
}

//----------------------------------------------------------------------------
vtkLabelInterpolator::~vtkLabelInterpolator()
{
  if (this->KernelLookupTable[0])
  {
    this->FreeKernelLookupTable();
  }
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "RadiusFactors: " << this->RadiusFactors[0] << " "
     << this->RadiusFactors[1] << " " << this->RadiusFactors[2] << "\n";
  os << indent << "BlurFactors: " << this->BlurFactors[0] << " "
     << this->BlurFactors[1] << " " << this->BlurFactors[2] << "\n";
  os << indent << "Antialiasing: "
     << (this->Antialiasing ? "On\n" : "Off\n");
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::ComputeSupportSize(
  const double matrix[16], int size[3])
{
  // compute the default support size for when matrix is null
  if (this->Antialiasing)
  {
    size[0] = VTK_LABEL_KERNEL_SIZE_MAX;
    size[1] = VTK_LABEL_KERNEL_SIZE_MAX;
    size[2] = VTK_LABEL_KERNEL_SIZE_MAX;
  }
  else
  {
    for (int i = 0; i < 3; i++)
    {
      // use blur factors to compute support size
      size[i] = 2*static_cast<int>(
        this->RadiusFactors[i] + 1.0 - VTK_INTERPOLATE_FLOOR_TOL);
      double rowscale = this->BlurFactors[i];
      if (rowscale > (1.0 + VTK_INTERPOLATE_FLOOR_TOL))
      {
        size[i] = 2*static_cast<int>(
          rowscale*this->RadiusFactors[i] + 1.0 - VTK_INTERPOLATE_FLOOR_TOL);
      }
    }
  }

  if (matrix == NULL)
  {
    return;
  }

  if (this->Antialiasing)
  {
    // if antialiasing is on, initialize blur factors to 1
    for (int i = 0; i < 3; i++)
    {
      this->BlurFactors[i] = 1.0;
      this->KernelSize[i] = 2*static_cast<int>(
        this->RadiusFactors[i] + 1.0 - VTK_INTERPOLATE_FLOOR_TOL);
    }
  }
  else
  {
    // keep blur factors, use kernel size computed from blur factors
    this->KernelSize[0] = size[0];
    this->KernelSize[1] = size[1];
    this->KernelSize[2] = size[2];
  }

  // if matrix does perspective, use the defaults just computed
  if (matrix[12] != 0 || matrix[13] != 0 || matrix[14] != 0 ||
      matrix[15] != 1.0)
  {
    return;
  }

  // use matrix to compute blur factors and kernel size
  for (int i = 0; i < 3; i++)
  {
    double rowscale = 0.0;
    for (int j = 0; j < 3; j++)
    {
      // compute the scale from a row of the matrix
      double x = matrix[4*i + j];
      rowscale += x*x;

      // verify that the element is an integer:
      // check fraction that remains after floor operation
      double f;
      vtkInterpolationMath::Floor(x, f);
    }

    if (this->Antialiasing)
    {
      // rowscale is the subsampling factor in a particular direction
      rowscale = sqrt(rowscale);
    }
    else
    {
      // ignore computed value, use factor provided by SetBlurFactors()
      rowscale = this->BlurFactors[i];
    }

    // if scale is greater than one, expand kernel size
    if (rowscale > (1.0 + VTK_INTERPOLATE_FLOOR_TOL))
    {
      // need extra suport for antialiasing
      this->BlurFactors[i] = rowscale;
      int s = 2*static_cast<int>(
        rowscale*this->RadiusFactors[i] + 1.0 - VTK_INTERPOLATE_FLOOR_TOL);
      size[i] = s;
      this->KernelSize[i] = s;
    }
  }

  // rebuild the kernel lookup tables
  this->InternalUpdate();
}

//----------------------------------------------------------------------------
bool vtkLabelInterpolator::IsSeparable()
{
  return true;
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::SetRadiusFactors(
  double x, double y, double z)
{
  if (this->RadiusFactors[0] != x ||
      this->RadiusFactors[1] != y ||
      this->RadiusFactors[2] != z)
  {
    this->RadiusFactors[0] = x;
    this->RadiusFactors[1] = y;
    this->RadiusFactors[2] = z;
    this->Modified();
  }
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::SetBlurFactors(double x, double y, double z)
{
  if (this->BlurFactors[0] != x ||
      this->BlurFactors[1] != y ||
      this->BlurFactors[2] != z)
  {
    this->BlurFactors[0] = x;
    this->BlurFactors[1] = y;
    this->BlurFactors[2] = z;
    this->Modified();
  }
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::SetAntialiasing(int val)
{
  val = (val != 0);
  if (this->Antialiasing != val)
  {
    this->Antialiasing = val;
    this->Modified();
  }
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::InternalDeepCopy(
  vtkAbstractImageInterpolator *a)
{
  vtkLabelInterpolator *obj = vtkLabelInterpolator::SafeDownCast(a);
  if (obj)
  {
    this->SetRadiusFactors(obj->RadiusFactors);
    this->SetAntialiasing(obj->Antialiasing);
    if (this->Antialiasing)
    {
      // if blur factors were computed, then don't call "modified"
      obj->GetBlurFactors(this->BlurFactors);
    }
    else
    {
      this->SetBlurFactors(obj->BlurFactors);
    }
  }

  this->KernelSize[0] = 6;
  this->KernelSize[1] = 6;
  this->KernelSize[2] = 6;

  if (this->KernelLookupTable[0])
  {
    this->FreeKernelLookupTable();
  }
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::InternalUpdate()
{
  bool blurchange = false;
  int mode = 0;
  int hsize[3];
  for (int i = 0; i < 3; i++)
  {
    static int minsize = 1;
    static int maxsize = VTK_LABEL_KERNEL_SIZE_MAX/2;
    int size = this->KernelSize[i]/2;
    size = ((size > minsize) ? size : minsize);
    size = ((size < maxsize) ? size : maxsize);
    hsize[i] = size;
    blurchange |= (fabs(this->BlurFactors[i] - this->LastBlurFactors[i]) >=
                   VTK_INTERPOLATE_FLOOR_TOL);
  }

  if (this->BlurFactors[0] > 1.0 + VTK_INTERPOLATE_FLOOR_TOL)
  {
    mode |= VTK_INTERPOLATION_WINDOW_XBLUR_MASK;
  }
  if (this->BlurFactors[1] > 1.0 + VTK_INTERPOLATE_FLOOR_TOL)
  {
    mode |= VTK_INTERPOLATION_WINDOW_YBLUR_MASK;
  }
  if (this->BlurFactors[2] > 1.0 + VTK_INTERPOLATE_FLOOR_TOL)
  {
    mode |= VTK_INTERPOLATION_WINDOW_ZBLUR_MASK;
  }

  mode |= (hsize[0] << VTK_INTERPOLATION_WINDOW_XSIZE_SHIFT);
  mode |= (hsize[1] << VTK_INTERPOLATION_WINDOW_YSIZE_SHIFT);
  mode |= (hsize[2] << VTK_INTERPOLATION_WINDOW_ZSIZE_SHIFT);

  if (this->InterpolationInfo->InterpolationMode != mode ||
      blurchange ||
      this->KernelLookupTable[0] == NULL)
  {
    this->BuildKernelLookupTable();
  }

  this->InterpolationInfo->InterpolationMode = mode;
  this->InterpolationInfo->ExtraInfo = this->KernelLookupTable;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//  Interpolation subroutines and associated code
//----------------------------------------------------------------------------

namespace {

//----------------------------------------------------------------------------
// Gaussian kernel computation: compute half of the interpolation
// kernel to make a lookup table of size "size".
// In the table, x=0.0 corresponds to index position zero, and
// x = 1.0 corresponds to index position "size", which is just
// beyond the end of the table and holds an implicit value of zero.

struct vtkGaussKernel
{
  template<class F>
  static void Evaluate(F *kernel, int size, double p);
};

template<class F>
void vtkGaussKernel::Evaluate(F *kernel, int size, double p)
{
  const double f = 2.5066282746310002; // sqrt(2*pi)
  double x = 0.0;
  do
  {
    *kernel++ = exp(-0.5*(x*x)*(f*f));
    x += p;
  }
  while (--size);
}

//----------------------------------------------------------------------------
template<class T, class F>
void vtkGaussInterpWeights(T *kernel, F *fX, F fx, int m)
{
  // table bins per unit
  int p = VTK_LABEL_KERNEL_TABLE_DIVISIONS;

  // compute table interpolation info
  F f = fx*p;
  int offset = static_cast<int>(f);
  f -= offset;
  F r = 1 - f;

  // interpolate the table
  int n = m;
  int i = (1 - (m >> 1))*p - offset;
  do
  {
    int i0 = i;
    int i1 = i + 1;
    int ni = -i0;
    i0 = ((i0 >= 0) ? i0 : ni);
    ni = -i1;
    i1 = ((i1 >= 0) ? i1 : ni);
    *fX++ = r*kernel[i0] + f*kernel[i1];
    i += p;
  }
  while (--n);
}

//----------------------------------------------------------------------------
template<class F, class T>
struct vtkImageLabelInterpolate
{
  static void General(
    vtkInterpolationInfo *info, const F point[3], F *outPtr);
};

//----------------------------------------------------------------------------
template <class F, class T>
void vtkImageLabelInterpolate<F, T>::General(
  vtkInterpolationInfo *info, const F point[3], F *outPtr)
{
  const T *inPtr = static_cast<const T *>(info->Pointer);
  int *inExt = info->Extent;
  vtkIdType *inInc = info->Increments;
  int numscalars = info->NumberOfComponents;

  // kernel lookup table
  float **kernel = static_cast<float **>(info->ExtraInfo);

  // size of kernel
  int mode = info->InterpolationMode;
  int xm = 2*((mode & VTK_INTERPOLATION_WINDOW_XSIZE_MASK)
              >> VTK_INTERPOLATION_WINDOW_XSIZE_SHIFT);
  int ym = 2*((mode & VTK_INTERPOLATION_WINDOW_YSIZE_MASK)
              >> VTK_INTERPOLATION_WINDOW_YSIZE_SHIFT);
  int zm = 2*((mode & VTK_INTERPOLATION_WINDOW_ZSIZE_MASK)
              >> VTK_INTERPOLATION_WINDOW_ZSIZE_SHIFT);

  // index to kernel midpoint position
  int xm2 = ((xm - 1) >> 1);
  int ym2 = ((ym - 1) >> 1);
  int zm2 = ((zm - 1) >> 1);

  F fx, fy, fz;
  int inIdX0 = vtkInterpolationMath::Floor(point[0], fx);
  int inIdY0 = vtkInterpolationMath::Floor(point[1], fy);
  int inIdZ0 = vtkInterpolationMath::Floor(point[2], fz);

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

  // the memory offsets
  vtkIdType factX[VTK_LABEL_KERNEL_SIZE_MAX];
  vtkIdType factY[VTK_LABEL_KERNEL_SIZE_MAX];
  vtkIdType factZ[VTK_LABEL_KERNEL_SIZE_MAX];

  // handle the borders
  int xi = inIdX0 - xm2;
  int yi = inIdY0 - ym2;
  int zi = inIdZ0 - zm2;
  int mm = xm;
  mm = ((mm >= ym) ? mm : ym);
  mm = ((mm >= zm) ? mm : zm);

  switch (info->BorderMode)
  {
    case VTK_IMAGE_BORDER_REPEAT:
    {
      int l = 0;
      do
      {
        factX[l] = vtkInterpolationMath::Wrap(xi, minX, maxX)*inIncX;
        factY[l] = vtkInterpolationMath::Wrap(yi, minY, maxY)*inIncY;
        factZ[l] = vtkInterpolationMath::Wrap(zi, minZ, maxZ)*inIncZ;
        l++; xi++; yi++; zi++;
      }
      while (--mm);
    }
      break;

    case VTK_IMAGE_BORDER_MIRROR:
    {
      int l = 0;
      do
      {
        factX[l] = vtkInterpolationMath::Mirror(xi, minX, maxX)*inIncX;
        factY[l] = vtkInterpolationMath::Mirror(yi, minY, maxY)*inIncY;
        factZ[l] = vtkInterpolationMath::Mirror(zi, minZ, maxZ)*inIncZ;
        l++; xi++; yi++; zi++;
      }
      while (--mm);
    }
      break;

    default:
    {
      int l = 0;
      do
      {
        factX[l] = vtkInterpolationMath::Clamp(xi, minX, maxX)*inIncX;
        factY[l] = vtkInterpolationMath::Clamp(yi, minY, maxY)*inIncY;
        factZ[l] = vtkInterpolationMath::Clamp(zi, minZ, maxZ)*inIncZ;
        l++; xi++; yi++; zi++;
      }
      while (--mm);
    }
      break;
  }

  // compute the kernel weights
  F fX[VTK_LABEL_KERNEL_SIZE_MAX];
  F fY[VTK_LABEL_KERNEL_SIZE_MAX];
  F fZ[VTK_LABEL_KERNEL_SIZE_MAX];

  vtkGaussInterpWeights(kernel[0], fX, fx, xm);
  vtkGaussInterpWeights(kernel[1], fY, fy, ym);
  vtkGaussInterpWeights(kernel[2], fZ, fz, zm);

  // the labels
  F labelWeights[VTK_LABEL_KERNEL_LABELS_MAX];
  T labels[VTK_LABEL_KERNEL_LABELS_MAX];

  // check if only one slice in a particular direction
  int multipleY = (minY != maxY);
  int multipleZ = (minZ != maxZ);

  // the limits to use when doing the interpolation
  int k1 = zm2*(1 - multipleZ);
  int k2 = (zm2 + 1)*(multipleZ + 1) - 1;
  int j1 = ym2*(1 - multipleY);
  int j2 = (ym2 + 1)*(multipleY + 1) - 1;

  do // loop over components
  {
    int labelCount = 0;
    int k = k1;
    do // loop over z
    {
      F ifz = fZ[k];
      vtkIdType factz = factZ[k];
      int j = j1;
      do // loop over y
      {
        F ify = fY[j];
        F fzy = ifz*ify;
        vtkIdType factzy = factz + factY[j];
        // loop over x
        const T *tmpPtr = inPtr + factzy;
        const F *tmpfX = fX;
        const vtkIdType *tmpfactX = factX;
        int l = xm;
        do
        {
          F val = fzy*(*tmpfX++);
          T label = tmpPtr[*tmpfactX++];
          int li = 0;
          for (; li < labelCount; li++)
          {
            if (label == labels[li])
            {
              labelWeights[li] += val;
              break;
            }
          }
          if (li == labelCount && li < VTK_LABEL_KERNEL_LABELS_MAX)
          {
            labels[li] = label;
            labelWeights[li] = val;
            labelCount++;
          }
        }
        while (--l);
      }
      while (++j <= j2);
    }
    while (++k <= k2);

    F maxweight = 0;
    T label = 0;
    for (int li = 0; li < labelCount; li++)
    {
      F weight = labelWeights[li];
      if (weight >= maxweight)
      {
        maxweight = weight;
        label = labels[li];
      }
    }
    *outPtr++ = label;
    inPtr++;
  }
  while (--numscalars);
}

//----------------------------------------------------------------------------
// Get the interpolation function for the specified data types
template<class F>
void vtkLabelInterpolatorGetInterpolationFunc(
  void (**interpolate)(vtkInterpolationInfo *, const F [3], F *),
  int dataType)
{
  switch (dataType)
  {
    vtkTemplateAliasMacro(
      *interpolate =
        &(vtkImageLabelInterpolate<F, VTK_TT>::General)
      );
    default:
      *interpolate = 0;
  }
}

//----------------------------------------------------------------------------
// Interpolation for precomputed weights

template <class F, class T>
struct vtkImageLabelRowInterpolate
{
  static void General(
    vtkInterpolationWeights *weights, int idX, int idY, int idZ,
    F *outPtr, int n);
};


//--------------------------------------------------------------------------
// helper function for high-order interpolation
template<class F, class T>
void vtkImageLabelRowInterpolate<F, T>::General(
  vtkInterpolationWeights *weights, int idX, int idY, int idZ,
  F *outPtr, int n)
{
  // the labels
  F labelWeights[VTK_LABEL_KERNEL_LABELS_MAX];
  T labels[VTK_LABEL_KERNEL_LABELS_MAX];

  int stepX = weights->KernelSize[0];
  int stepY = weights->KernelSize[1];
  int stepZ = weights->KernelSize[2];
  idX *= stepX;
  idY *= stepY;
  idZ *= stepZ;
  const F *fX = static_cast<F *>(weights->Weights[0]) + idX;
  const F *fY = static_cast<F *>(weights->Weights[1]) + idY;
  const F *fZ = static_cast<F *>(weights->Weights[2]) + idZ;
  const vtkIdType *factX = weights->Positions[0] + idX;
  const vtkIdType *factY = weights->Positions[1] + idY;
  const vtkIdType *factZ = weights->Positions[2] + idZ;
  const T *inPtr = static_cast<const T *>(weights->Pointer);

  int numscalars = weights->NumberOfComponents;
  for (int i = n; i > 0; --i)
  {
    const T *inPtr0 = inPtr;
    int c = numscalars;
    do // loop over components
    {
      int labelCount = 0;
      int k = 0;
      do // loop over z
      {
        F ifz = fZ[k];
        vtkIdType factz = factZ[k];
        int j = 0;
        do // loop over y
        {
          F ify = fY[j];
          F fzy = ifz*ify;
          vtkIdType factzy = factz + factY[j];
          // loop over x
          const T *tmpPtr = inPtr0 + factzy;
          const F *tmpfX = fX;
          const vtkIdType *tmpfactX = factX;
          int l = stepX;
          do
          {
            F val = fzy*(*tmpfX++);
            T label = tmpPtr[*tmpfactX++];
            int li = 0;
            for (; li < labelCount; li++)
            {
              if (label == labels[li])
              {
                labelWeights[li] += val;
                break;
              }
            }
            if (li == labelCount && li < VTK_LABEL_KERNEL_LABELS_MAX)
            {
              labels[li] = label;
              labelWeights[li] = val;
              labelCount++;
            }
          }
          while (--l);
        }
        while (++j < stepY);
      }
      while (++k < stepZ);

      F maxweight = 0;
      T label = 0;
      for (int li = 0; li < labelCount; li++)
      {
        F weight = labelWeights[li];
        if (weight >= maxweight)
        {
          maxweight = weight;
          label = labels[li];
        }
      }
      *outPtr++ = label;
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
void vtkLabelInterpolatorGetRowInterpolationFunc(
  void (**summation)(vtkInterpolationWeights *weights, int idX, int idY,
                     int idZ, F *outPtr, int n),
  int scalarType)
{
  switch (scalarType)
  {
    vtkTemplateAliasMacro(
      *summation = &(vtkImageLabelRowInterpolate<F,VTK_TT>::General)
      );
    default:
      *summation = 0;
  }
}

//----------------------------------------------------------------------------
template<class F>
void vtkLabelInterpolatorPrecomputeWeights(
  const F newmat[16], const int outExt[6], int clipExt[6],
  const F bounds[6], vtkInterpolationWeights *weights)
{
  float **kernel = static_cast<float **>(weights->ExtraInfo);
  weights->WeightType = vtkTypeTraits<F>::VTKTypeID();
  int sizes[3];
  bool blur[3];
  int mode = weights->InterpolationMode;
  sizes[0] = 2*((mode & VTK_INTERPOLATION_WINDOW_XSIZE_MASK)
                >> VTK_INTERPOLATION_WINDOW_XSIZE_SHIFT);
  sizes[1] = 2*((mode & VTK_INTERPOLATION_WINDOW_YSIZE_MASK)
                >> VTK_INTERPOLATION_WINDOW_YSIZE_SHIFT);
  sizes[2] = 2*((mode & VTK_INTERPOLATION_WINDOW_ZSIZE_MASK)
                >> VTK_INTERPOLATION_WINDOW_ZSIZE_SHIFT);
  blur[0] = ((mode & VTK_INTERPOLATION_WINDOW_XBLUR_MASK) != 0);
  blur[1] = ((mode & VTK_INTERPOLATION_WINDOW_YBLUR_MASK) != 0);
  blur[2] = ((mode & VTK_INTERPOLATION_WINDOW_ZBLUR_MASK) != 0);

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

    // the kernel size should not exceed the input dimension
    int m = sizes[j];
    int m2 = ((m - 1) >> 1);
    int step = m;
    int inCount = maxExt - minExt + 1;
    step = ((step < inCount) ? step : inCount);

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

      F f = 0;
      int idx = vtkInterpolationMath::Floor(point, f);
      int lmax = 1;
      if (step > 1)
      {
        idx -= m2;
        lmax = m;
      }

      int inId[VTK_LABEL_KERNEL_SIZE_MAX];

      int l = 0;
      switch (weights->BorderMode)
      {
        case VTK_IMAGE_BORDER_REPEAT:
          do
          {
            inId[l] = vtkInterpolationMath::Wrap(idx++, minExt, maxExt);
          }
          while (++l < lmax);
          break;

        case VTK_IMAGE_BORDER_MIRROR:
          do
          {
            inId[l] = vtkInterpolationMath::Mirror(idx++, minExt, maxExt);
          }
          while (++l < lmax);
          break;

        default:
           do
           {
            inId[l] = vtkInterpolationMath::Clamp(idx++, minExt, maxExt);
           }
          while (++l < lmax);
          break;
      }

      // compute the weights and offsets
      vtkIdType inInc = weights->Increments[k];
      if (step == 1)
      {
        positions[step*i] = inId[0]*inInc;
        constants[step*i] = static_cast<F>(1);
      }
      else
      {
        F g[VTK_LABEL_KERNEL_SIZE_MAX];
        vtkGaussInterpWeights(kernel[j], g, f, m);

        if (step == m)
        {
          int ll = 0;
          do
          {
            positions[step*i + ll] = inId[ll]*inInc;
            constants[step*i + ll] = g[ll];
          }
          while (++ll < step);
        }
        else
        {
          // it gets tricky if the data is thinner than the kernel
          F gg[VTK_LABEL_KERNEL_SIZE_MAX];
          int ll = 0;
          do { gg[ll] = 0; } while (++ll < m);
          ll = 0;
          do
          {
            int rIdx = inId[ll];
            gg[rIdx] += g[ll];
          }
          while (++ll < m);
          ll = 0;
          do
          {
            positions[step*i + ll] = ll*inInc;
            constants[step*i + ll] = gg[ll];
          }
          while (++ll < step);
        }
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
void vtkLabelInterpolator::GetInterpolationFunc(
  void (**func)(vtkInterpolationInfo *, const double [3], double *))
{
  vtkLabelInterpolatorGetInterpolationFunc(
    func, this->InterpolationInfo->ScalarType);
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::GetInterpolationFunc(
  void (**func)(vtkInterpolationInfo *, const float [3], float *))
{
  vtkLabelInterpolatorGetInterpolationFunc(
    func, this->InterpolationInfo->ScalarType);
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::GetRowInterpolationFunc(
  void (**func)(vtkInterpolationWeights *, int, int, int, double *, int))
{
  vtkLabelInterpolatorGetRowInterpolationFunc(
    func, this->InterpolationInfo->ScalarType);
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::GetRowInterpolationFunc(
  void (**func)(vtkInterpolationWeights *, int, int, int, float *, int))
{
  vtkLabelInterpolatorGetRowInterpolationFunc(
    func, this->InterpolationInfo->ScalarType);
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::PrecomputeWeightsForExtent(
  const double matrix[16], const int extent[6], int newExtent[6],
  vtkInterpolationWeights *&weights)
{
  weights = new vtkInterpolationWeights(*this->InterpolationInfo);

  vtkLabelInterpolatorPrecomputeWeights(
    matrix, extent, newExtent, this->StructuredBoundsDouble, weights);
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::PrecomputeWeightsForExtent(
  const float matrix[16], const int extent[6], int newExtent[6],
  vtkInterpolationWeights *&weights)
{
  weights = new vtkInterpolationWeights(*this->InterpolationInfo);

  vtkLabelInterpolatorPrecomputeWeights(
    matrix, extent, newExtent, this->StructuredBoundsFloat, weights);
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::FreePrecomputedWeights(
  vtkInterpolationWeights *&weights)
{
  this->Superclass::FreePrecomputedWeights(weights);
}

//----------------------------------------------------------------------------
// build any tables required for the interpolation
void vtkLabelInterpolator::BuildKernelLookupTable()
{
  if (this->KernelLookupTable[0])
  {
    this->FreeKernelLookupTable();
  }

  float *kernel[3];
  kernel[0] = 0;
  kernel[1] = 0;
  kernel[2] = 0;

  for (int i = 0; i < 3; i++)
  {
    // reuse the X kernel lookup table if possible
    if (i > 0 && this->KernelSize[i] == this->KernelSize[0] &&
        fabs(this->RadiusFactors[i] - this->RadiusFactors[0]) <
          VTK_INTERPOLATE_FLOOR_TOL &&
        fabs(this->BlurFactors[i] - this->BlurFactors[0]) <
          VTK_INTERPOLATE_FLOOR_TOL)
    {
      kernel[i] = kernel[0];
      continue;
    }

    // kernel parameters
    int m = this->KernelSize[i];
    double b = this->BlurFactors[i];

    // blur factor must be restricted to half the max kernel size
    if (b > 0.5*VTK_LABEL_KERNEL_SIZE_MAX)
    {
      b = 0.5*VTK_LABEL_KERNEL_SIZE_MAX;
    }

    // compute lookup table size and step size
    int size = m/2*VTK_LABEL_KERNEL_TABLE_DIVISIONS;
    double p = 1.0/(b*VTK_LABEL_KERNEL_TABLE_DIVISIONS);

    // allocate and compute the kernel lookup table
    // (add a small safety buffer that will be filled with zeros)
    kernel[i] = new float[size + 4];

    int cutoff = static_cast<int>(this->RadiusFactors[i]*b/p + 0.5);
    cutoff = (cutoff < size ? cutoff : size);
    cutoff += 1;
    vtkGaussKernel::Evaluate(kernel[i], cutoff, p);

    // add a tail of zeros for when table is interpolated
    float *kptr = &kernel[i][cutoff];
    int k = size + 4 - cutoff;
    do
    {
      *kptr++ = 0;
    }
    while (--k);

    if (b > 1.0)
    {
      // if kernel stretched to create blur, divide by stretch factor
      float *ktmp = kernel[i];
      float bf = 1.0/b;
      int j = size;
      do
      {
        *ktmp *= bf;
        ktmp++;
      }
      while (--j);
    }
  }

  this->KernelLookupTable[0] = kernel[0];
  this->KernelLookupTable[1] = kernel[1];
  this->KernelLookupTable[2] = kernel[2];

  this->LastBlurFactors[0] = this->BlurFactors[0];
  this->LastBlurFactors[1] = this->BlurFactors[1];
  this->LastBlurFactors[2] = this->BlurFactors[2];
}

//----------------------------------------------------------------------------
void vtkLabelInterpolator::FreeKernelLookupTable()
{
  float *kernel = this->KernelLookupTable[0];
  if (kernel)
  {
    delete [] kernel;
    for (int i = 1; i < 3; i++)
    {
      if (this->KernelLookupTable[i] != kernel)
      {
        delete [] this->KernelLookupTable[i];
      }
    }
  }
}
