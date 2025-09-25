/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImageAmoebaGrid.cxx,v $
  Language:  C++
  Date:      $Date: 2007/08/24 20:02:25 $
  Version:   $Revision: 1.8 $

Copyright (c) 1993-2000 Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
   of any contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 * Modified source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#include "vtkImageAmoebaGrid.h"

#include "vtkImageData.h"
#include "vtkImageStencilData.h"
#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkExecutive.h"

//----------------------------------------------------------------------------
vtkImageAmoebaGrid* vtkImageAmoebaGrid::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkImageAmoebaGrid");
  if(ret)
  {
    return (vtkImageAmoebaGrid*)ret;
  }
  // If the factory was unable to create the object, then create it here.
  return new vtkImageAmoebaGrid;
}

// on i386 platforms, SetOptimization(2) provides integer-based calculations
// (on other platforms, integer math is slower than float math)
#ifdef i386  // kwang to test the performance gain on altix IA64 system
#define VTK_RESLICE_INTEGER_MATH 1
#endif

//----------------------------------------------------------------------------
vtkImageAmoebaGrid::vtkImageAmoebaGrid()
{
  this->ShrinkFactors[0] = 1.0f;
  this->ShrinkFactors[1] = 1.0f;
  this->ShrinkFactors[2] = 1.0f;

  this->KernelRadius[0] = 1;
  this->KernelRadius[1] = 1;
  this->KernelRadius[2] = 1;

  this->ReverseStencil = 0;

  this->LastThreadCount = 0;
  this->VectorLength = NULL;
  this->VectorsMinimized = NULL;
  this->TotalCost = NULL;

  this->Tolerance = 0.005;

  // we have the image inputs and the optional stencil input
  this->SetNumberOfInputPorts(2);
}

//----------------------------------------------------------------------------
vtkImageAmoebaGrid::~vtkImageAmoebaGrid()
{
  if (this->VectorLength)
  {
    delete [] this->VectorLength;
  }

  if (this->VectorsMinimized)
  {
    delete [] this->VectorsMinimized;
  }

  if (this->TotalCost)
  {
    delete [] this->TotalCost;
  }
}

//----------------------------------------------------------------------------
int vtkImageAmoebaGrid::FillInputPortInformation(int port,
                                                 vtkInformation *info)
{
  if (port == 1)
  {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageStencilData");
    // the stencil input is optional
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  }
  else
  {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
  }
  return 1;
}

//----------------------------------------------------------------------------
void vtkImageAmoebaGrid::SetStencil(vtkImageStencilData *stencil)
{
  this->SetNthInputConnection(1, 0,
    (stencil ? stencil->GetProducerPort() : 0));
}

//----------------------------------------------------------------------------
vtkImageStencilData *vtkImageAmoebaGrid::GetStencil()
{
  if (this->GetNumberOfInputConnections(1) < 1)
  {
    return NULL;
  }
  return vtkImageStencilData::SafeDownCast(
    this->GetExecutive()->GetInputData(1, 0));
}

#ifdef VTK_RESLICE_INTEGER_MATH
//----------------------------------------------------------------------------
// for fixed-point math, define a basic set of macros

// the radix is 14, to support signed values up to 1<<16
#define VTK_FP_RADIX 14
#define VTK_FP_RADIX_MINUS_1 13

// 0.5 in fixed-point
#define VTK_FP_HALF (1<<VTK_FP_RADIX_MINUS_1)

// various integer values in fixed-point
#define VTK_FP_0 0
#define VTK_FP_1 (1<<VTK_FP_RADIX)
#define VTK_FP_2 (2<<VTK_FP_RADIX)
#define VTK_FP_3 (3<<VTK_FP_RADIX)
#define VTK_FP_4 (4<<VTK_FP_RADIX)

// a very nifty hack I discovered at http://www.stereopsis.com/FPU.html,
// it adds (2**(52-radix))*1.5 to a number to get just the right roundoff
// from double to fixed-point
static inline int vtkCastFloatToFixed(double x)
{
  union { double d; unsigned int i[2]; } dual;
  dual.d = x + 412316860416.0; // (2**(52-radix))*1.5
#ifdef VTK_WORDS_BIGENDIAN
  return dual.i[1];
#else
  return dual.i[0];
#endif
}

//----------------------------------------------------------------------------
// converting the other way is much more straightforward
static inline double vtkCastFixedToFloat(int x)
{
  return x*(1.0/VTK_FP_1);
}

//----------------------------------------------------------------------------
// what follows are a whole bunch of inline functions that provide
// equivalent behaviour for fixed-point numbers vs. floats (note
// that doubles are not supported, but would be easy to add)
static inline int vtkResliceFloor(int x)
{
  return x>>VTK_FP_RADIX;
}

//----------------------------------------------------------------------------
static inline int vtkResliceCeil(int x)
{
  return ((1<<VTK_FP_RADIX) - 1 + x)>>VTK_FP_RADIX;
}

//----------------------------------------------------------------------------
static inline int vtkResliceRound(int x)
{
  return (x + VTK_FP_HALF)>>VTK_FP_RADIX;
}

//----------------------------------------------------------------------------
// convert a fixed-point into an integer plus a fraction
static inline int vtkResliceFloor(int x, int &f)
{
  int ix = x>>VTK_FP_RADIX;
  f = x - (ix<<VTK_FP_RADIX);

  return ix;
}

//----------------------------------------------------------------------------
// multiplication, the product must be less than 8 for fixed-point
static inline int vtkResliceQuikMul(int xy)
{
  return (xy + VTK_FP_HALF)>>VTK_FP_RADIX;
}

//----------------------------------------------------------------------------
// multiplication of larger fixed-point numbers, does not check overflow
static inline int vtkResliceMultiply(int x, int y)
{
  int hx = x>>VTK_FP_RADIX;
  int hy = y>>VTK_FP_RADIX;
  int lx = x - (hx<<VTK_FP_RADIX);
  int ly = y - (hy<<VTK_FP_RADIX);

  return ((lx*ly + VTK_FP_HALF)>>VTK_FP_RADIX) + hx*ly + x*hy;
}

//----------------------------------------------------------------------------
static inline int vtkResliceInverse(int x)
{
  return ((1<<(2*VTK_FP_RADIX + 1))/x + 1)>>1;
}

//----------------------------------------------------------------------------
// cast between float and fixed-point
static inline void vtkResliceCast(float x, int& y)
{
  y = vtkCastFloatToFixed(x);
}

//----------------------------------------------------------------------------
static inline void vtkResliceCast(double x, int& y)
{
  y = vtkCastFloatToFixed(x);
}

//----------------------------------------------------------------------------
static inline void vtkResliceCast(int x, float& y)
{
  y = vtkCastFixedToFloat(x);
}

//----------------------------------------------------------------------------
// 1-x, it gets used a lot
static inline int vtkResliceOneMinusX(int x)
{
  return VTK_FP_1 - x;
}

//----------------------------------------------------------------------------
// check if a number is equal to one
static inline int vtkResliceIsEqualToOne(int x)
{
  return (VTK_FP_1 == x);
}

//----------------------------------------------------------------------------
// check if a number is an integer
static inline int vtkResliceIsInteger(int x)
{
  return (x == ((x>>VTK_FP_RADIX)<<VTK_FP_RADIX));
}

//----------------------------------------------------------------------------
// rounding functions for fixed-point, note that only
// char, unsigned char, short, unsigned short are supported
static inline void vtkResliceRound(int val, char& rnd)
{
  rnd = (char)vtkResliceRound(val);
}

//----------------------------------------------------------------------------
static inline void vtkResliceRound(int val, unsigned char& rnd)
{
  rnd = (unsigned char)vtkResliceRound(val);
}

//----------------------------------------------------------------------------
static inline void vtkResliceRound(int val, short& rnd)
{
  rnd = (short)vtkResliceRound(val);
}

//----------------------------------------------------------------------------
static inline void vtkResliceRound(int val, unsigned short& rnd)
{
  rnd = (unsigned short)vtkResliceRound(val);
}
#endif /* VTK_RESLICE_INTEGER_MATH */

//----------------------------------------------------------------------------
// fast floor() function for converting a float to an int
// (the floor() implementation on some computers is much slower than this,
// because they require some 'exact' behaviour that we don't).

// The 'floor' function on x86 and mips is many times slower than these
// and is used a lot in this code, optimize for different CPU architectures
inline int vtkResliceFloor(double x)
{
#if defined mips || defined sparc || defined __ppc__
  x += 2147483648.0;
  unsigned int i = (unsigned int)(x);
  return (int)(i - 2147483648U);
#elif defined i386 || defined _M_IX86
  union { double d; unsigned short s[4]; unsigned int i[2]; } dual;
  dual.d = x + 103079215104.0;  // (2**(52-16))*1.5
  return (int)((dual.i[1]<<16)|((dual.i[0])>>16));
#elif defined ia64 || defined __ia64__ || defined IA64
  x += 103079215104.0;
  long long i = (long long)(x);
  return (int)(i - 103079215104LL);
#else
  double y = floor(x);
  return (int)(y);
#endif
}

//----------------------------------------------------------------------------
static inline int vtkResliceCeil(double x)
{
  return -vtkResliceFloor(-x - 1.0) - 1;
}

//----------------------------------------------------------------------------
static inline int vtkResliceRound(double x)
{
  return vtkResliceFloor(x + 0.5);
}

//----------------------------------------------------------------------------
static inline int vtkResliceFloor(float x)
{
  return vtkResliceFloor((double)x);
}

//----------------------------------------------------------------------------
static inline int vtkResliceCeil(float x)
{
  return vtkResliceCeil((double)x);
}

//----------------------------------------------------------------------------
static inline int vtkResliceRound(float x)
{
  return vtkResliceRound((double)x);
}

//----------------------------------------------------------------------------
// convert a float into an integer plus a fraction
static inline int vtkResliceFloor(float x, float &f)
{
  int ix = vtkResliceFloor(x);
  f = x - ix;

  return ix;
}

//----------------------------------------------------------------------------
static inline float vtkResliceQuikMul(float xy)
{
  return xy;
}

//----------------------------------------------------------------------------
static inline float vtkResliceMultiply(float x, float y)
{
  return x*y;
}

//----------------------------------------------------------------------------
// invert a number
static inline float vtkResliceInverse(float x)
{
  return 1.0f/x;
}

//----------------------------------------------------------------------------
static inline void vtkResliceCast(float x, float& y)
{
  y = x;
}

//----------------------------------------------------------------------------
static inline void vtkResliceCast(double x, float& y)
{
  y = (float)x;
}

//----------------------------------------------------------------------------
static inline float vtkResliceOneMinusX(float x)
{
  return 1.0f - x;
}

//----------------------------------------------------------------------------
static inline int vtkResliceIsEqualToOne(float x)
{
  return (1.0f == x);
}

//----------------------------------------------------------------------------
static inline int vtkResliceIsInteger(float x)
{
  return (x == vtkResliceFloor(x));
}


//----------------------------------------------------------------------------
// rounding functions for each type, with some crazy stunts to avoid
// the use of the 'floor' function which is too slow on x86
template<class T>
static inline void vtkResliceRound(float val, T& rnd)
{
  rnd = vtkResliceRound(val);
}

//----------------------------------------------------------------------------
template<class T>
static inline void vtkResliceRound(double val, T& rnd)
{
  rnd = vtkResliceRound(val);
}

//----------------------------------------------------------------------------
static inline void vtkResliceRound(float val, float& rnd)
{
  rnd = val;
}

//----------------------------------------------------------------------------
static inline void vtkResliceRound(double val, float& rnd)
{
  rnd = val;
}

//----------------------------------------------------------------------------
static inline void vtkResliceRound(float val, double& rnd)
{
  rnd = val;
}

//----------------------------------------------------------------------------
static inline void vtkResliceRound(double val, double& rnd)
{
  rnd = val;
}

//----------------------------------------------------------------------------
void vtkImageAmoebaGrid::ComputeInputUpdateExtents(vtkDataObject *output)
{
  this->vtkImageMultipleInputFilter::ComputeInputUpdateExtents(output);

  vtkImageStencilData *stencil = this->GetStencil();
  if (stencil)
  {
    stencil->SetUpdateExtent(output->GetUpdateExtent());
  }

  if (this->LastThreadCount != this->GetNumberOfThreads())
  {
    this->LastThreadCount = this->GetNumberOfThreads();
    if (this->VectorLength)
    {
      delete [] this->VectorLength;
    }
    this->VectorLength = new double[this->LastThreadCount];

    if (this->VectorsMinimized)
    {
      delete [] this->VectorsMinimized;
    }
    this->VectorsMinimized = new int[this->LastThreadCount];

    if (this->TotalCost)
    {
      delete [] this->TotalCost;
    }
    this->TotalCost = new double[this->LastThreadCount];

    // these arrays should be initialized otherwise will cause problem
    // during summing. add by kwang 05/20/2004
    for (int i=0; i<this->LastThreadCount; i++)
    {
        this->VectorLength[i] = 0.0f;
        this->VectorsMinimized[i] = 0;
        this->TotalCost[i] = 0.0f;
    }
  }
}

//----------------------------------------------------------------------------
// Compute and return the mean vector length
float vtkImageAmoebaGrid::GetMeanVectorLength()
{
  double totalLength = 0;
  for (int x=0; x<this->LastThreadCount; x++)
  {
    totalLength += this->VectorLength[x];
  }
  int totalVectors = this->GetVectorsMinimized();
  if (totalVectors)
  {
    return totalLength / totalVectors;
  }
  else
  {
    return 0.0f;
  }
}

//----------------------------------------------------------------------------
// Compute and return the mean cost function
float vtkImageAmoebaGrid::GetMeanCost()
{
  double totalCost = 0.0f;
  for (int x=0; x<this->LastThreadCount; x++)
  {
    totalCost += this->TotalCost[x];
  }
  int totalVectors = this->GetVectorsMinimized();
  if (totalVectors)
  {
    return totalCost / totalVectors;
  }
  else
  {
    return 0.0f;
  }
}

//----------------------------------------------------------------------------
// Compute and return the total number of vectors minimized in the most recent
// iteration.
int vtkImageAmoebaGrid::GetVectorsMinimized()
{
  int minimized = 0;
  for (int x=0; x<this->LastThreadCount; x++)
  {
    minimized += this->VectorsMinimized[x];
  }
  return minimized;
}

//----------------------------------------------------------------------------
// This method computes the Region of input necessary to generate outRegion.
// Computes a different region for the images than for the input grid since
// they could have different extents. We grab all of the input images. This
// could be made smarter, but for now we always want it all anyway.
void vtkImageAmoebaGrid::ComputeInputUpdateExtent(int inExt[6],
						  int outExt[6],
						  int whichInput)
{
  int *wholeInExt = this->GetInput(whichInput)->GetWholeExtent();

  if (whichInput>2)//The input grid or the Stencil
  {
    memcpy(inExt, outExt, sizeof(int)*6);
  }
  else // one of the two images
  {
    memcpy(inExt,wholeInExt, sizeof(int)*6);
  }
}

//----------------------------------------------------------------------------
//Set Up the output to match the third input (0,1,2) extent
// inDatas[0] = Patient
// inDatas[1] = Model
// inDatas[2] = Input Grid
// Force 3 component floats for the output grid
void vtkImageAmoebaGrid::ExecuteInformation(vtkImageData **inDatas,
					    vtkImageData *outData)
{
  int wholeDataExt[6];
  int wholeGridExt[6];
  vtkFloatingPointType inDataSpacing[3];

  inDatas[0]->UpdateInformation();
  inDatas[1]->UpdateInformation();
  inDatas[2]->UpdateInformation();

  inDatas[1]->GetSpacing(inDataSpacing);
  inDatas[1]->GetWholeExtent(wholeDataExt);
  inDatas[2]->GetWholeExtent(wholeGridExt);

  this->ShrinkFactors[0] = ((float)wholeDataExt[1] - (float)wholeDataExt[0])/
                           ((float)wholeGridExt[1] - (float)wholeGridExt[0]);
  this->ShrinkFactors[1] = ((float)wholeDataExt[3] - (float)wholeDataExt[2])/
                           ((float)wholeGridExt[3] - (float)wholeGridExt[2]);
  this->ShrinkFactors[2] = ((float)wholeDataExt[5] - (float)wholeDataExt[4])/
                           ((float)wholeGridExt[5] - (float)wholeGridExt[4]);

  inDataSpacing[0]*=this->ShrinkFactors[0];
  inDataSpacing[1]*=this->ShrinkFactors[1];
  inDataSpacing[2]*=this->ShrinkFactors[2];

  outData->SetWholeExtent(wholeGridExt);
  outData->SetSpacing(inDataSpacing);
  outData->SetOrigin(inDatas[0]->GetOrigin());
  outData->SetNumberOfScalarComponents(3);
  outData->SetScalarType(VTK_FLOAT);

  // need to set the spacing and origin of the stencil to match the output
  vtkImageStencilData *stencil = this->GetStencil();
  if (stencil)
  {
    stencil->SetSpacing(inDataSpacing);
    stencil->SetOrigin(inDatas[0]->GetOrigin());
  }
}

//----------------------------------------------------------------------------
class _vtkAmoebaParms {
public:
  void *modelPtr;
  float modelExtent[6];
  vtkIdType modIncs[3];
  void *patPtr;
  int *patWholeInExt;
  vtkIdType patIncs[3];
  int patext[3]; // precompute patWholeInExt([1]-[0],[3]-[2],[5]-[4])
  int scalarType;
  float *hint;
  vtkFloatingPointType *spacing;
  int maxLength;
  float kernelRadius[3];
  int kernelDiameter;
  vtkFunctionMinimizer *Minimizer;
  double *disps;  // pointer to scalar vars used by Minimizer
  vtkImageData *patientData;
  vtkImageData *modelData;
  unsigned char *modelPoints;
  float modelVector[3]; // location of current vector in model data coordinates
  double sqrt_mod_sum_squared;
};

//----------------------------------------------------------------------------
template <class T>
static inline int GetModelBlock(T *modelPtr, _vtkAmoebaParms *pb)
{
  int idxX, idxY, idxZ, kernelDiameter;

  kernelDiameter = idxX = idxY = idxZ = (int)((pb->kernelRadius[0]+0.5)*2.0);

  pb->modelExtent[0] = pb->modelVector[0]-pb->kernelRadius[0];
  pb->modelExtent[1] = pb->modelVector[0]+pb->kernelRadius[0];
  pb->modelExtent[2] = pb->modelVector[1]-pb->kernelRadius[1];
  pb->modelExtent[3] = pb->modelVector[1]+pb->kernelRadius[1];
  pb->modelExtent[4] = pb->modelVector[2]-pb->kernelRadius[2];
  pb->modelExtent[5] = pb->modelVector[2]+pb->kernelRadius[2];

  int *modelWholeExt = pb->modelData->GetWholeExtent();
  vtkIdType *modelIncs = pb->modIncs;

  float fx,fy,fz;
  int floorX = vtkResliceFloor(pb->modelExtent[0],fx);
  int floorY = vtkResliceFloor(pb->modelExtent[2],fy);
  int floorZ = vtkResliceFloor(pb->modelExtent[4],fz);

  int gridId0X = floorX - modelWholeExt[0];
  int gridId0Y = floorY - modelWholeExt[2];
  int gridId0Z = floorZ - modelWholeExt[4];

  int gridId1X = gridId0X + 1;
  int gridId1Y = gridId0Y + 1;
  int gridId1Z = gridId0Z + 1;

  int extX = modelWholeExt[1]-modelWholeExt[0];
  int extY = modelWholeExt[3]-modelWholeExt[2];
  int extZ = modelWholeExt[5]-modelWholeExt[4];

  // do bounds check, most points will be inside so optimize for that
  if ((gridId0X | (extX - (gridId1X+kernelDiameter-1)) |
       gridId0Y | (extY - (gridId1Y+kernelDiameter-1)) |
       gridId0Z | (extZ - (gridId1Z+kernelDiameter-1))) < 0)
  {
    return 0;
  }
  // do trilinear interpolation
  vtkIdType factX0 = gridId0X*modelIncs[0];
  vtkIdType factY0 = gridId0Y*modelIncs[1];
  vtkIdType factZ0 = gridId0Z*modelIncs[2];

  vtkIdType factX1 = gridId1X*modelIncs[0];
  vtkIdType factY1 = gridId1Y*modelIncs[1];
  vtkIdType factZ1 = gridId1Z*modelIncs[2];

  T *p000 = modelPtr+factX0+factY0+factZ0;
  T *p001 = modelPtr+factX0+factY0+factZ1;
  T *p010 = modelPtr+factX0+factY1+factZ0;
  T *p011 = modelPtr+factX0+factY1+factZ1;
  T *p100 = modelPtr+factX1+factY0+factZ0;
  T *p101 = modelPtr+factX1+factY0+factZ1;
  T *p110 = modelPtr+factX1+factY1+factZ0;
  T *p111 = modelPtr+factX1+factY1+factZ1;

  float rx = 1.0 - fx;
  float ry = 1.0 - fy;
  float rz = 1.0 - fz;

  float ryrz = ry*rz;
  float ryfz = ry*fz;
  float fyrz = fy*rz;
  float fyfz = fy*fz;

  unsigned int irxryrz = (unsigned int)(rx*ryrz *512);
  unsigned int irxryfz = (unsigned int)(rx*ryfz *512);
  unsigned int irxfyrz = (unsigned int)(rx*fyrz *512);
  unsigned int irxfyfz = (unsigned int)(rx*fyfz *512);
  unsigned int ifxryrz = (unsigned int)(fx*ryrz *512);
  unsigned int ifxryfz = (unsigned int)(fx*ryfz *512);
  unsigned int ifxfyrz = (unsigned int)(fx*fyrz *512);
  unsigned int ifxfyfz = (unsigned int)(fx*fyfz *512);

  vtkIdType modContIncY = modelIncs[1] - (kernelDiameter);
  vtkIdType modContIncZ = modelIncs[2] - (kernelDiameter)*modelIncs[1];
  long int modData;

  // the array holding all the interpolated values
  unsigned char *modelPoints = pb->modelPoints;
  double value1, value2;
  pb->sqrt_mod_sum_squared = 0.0;
  do
  {
    do
    {
      do
      {
	modData = (irxryrz* *p000++ + irxryfz* *p001++ +
		   irxfyrz* *p010++ + irxfyfz* *p011++ +
		   ifxryrz* *p100++ + ifxryfz* *p101++ +
		   ifxfyrz* *p110++ + ifxfyfz* *p111++);
	*modelPoints = modData >>=9; // divide by 512
	//pb->sqrt_mod_sum_squared += *modelPoints * *modelPoints++;
	value1 = (double)*modelPoints;
	value2 = (double)*modelPoints++;
	pb->sqrt_mod_sum_squared += value1*value2;
      }
      while (--idxX);
      p000 += modContIncY; p001 += modContIncY;
      p010 += modContIncY; p011 += modContIncY;
      p100 += modContIncY; p101 += modContIncY;
      p110 += modContIncY; p111 += modContIncY;
      idxX = kernelDiameter;
    }
    while (--idxY);
    p000 += modContIncZ; p001 += modContIncZ;
    p010 += modContIncZ; p011 += modContIncZ;
    p100 += modContIncZ; p101 += modContIncZ;
    p110 += modContIncZ; p111 += modContIncZ;
    idxY = kernelDiameter;
  }
  while (--idxZ);
  pb->sqrt_mod_sum_squared = sqrt(pb->sqrt_mod_sum_squared);
  return 1;
}

//----------------------------------------------------------------------------
template <class T>
static inline T GetDataAtPoint(T *dataPtr, vtkIdType *dataIncs,
			       int *dataWholeExt, float *xyz)
{
  // change point into integer plus fraction
  float fx,fy,fz;
  int floorX = vtkResliceFloor(xyz[0],fx);
  int floorY = vtkResliceFloor(xyz[1],fy);
  int floorZ = vtkResliceFloor(xyz[2],fz);

  int gridId0X = floorX - dataWholeExt[0];
  int gridId0Y = floorY - dataWholeExt[2];
  int gridId0Z = floorZ - dataWholeExt[4];

  int gridId1X = gridId0X + 1;
  int gridId1Y = gridId0Y + 1;
  int gridId1Z = gridId0Z + 1;

  int extX = dataWholeExt[1]-dataWholeExt[0];
  int extY = dataWholeExt[3]-dataWholeExt[2];
  int extZ = dataWholeExt[5]-dataWholeExt[4];

  // do bounds check, most points will be inside so optimize for that
  if ((gridId0X | (extX - gridId1X) |
       gridId0Y | (extY - gridId1Y) |
       gridId0Z | (extZ - gridId1Z)) < 0)
  {
    if (gridId0X < 0)
    {
      gridId0X = 0;
      gridId1X = 0;
      fx = 0;
    }
    else if (gridId1X > extX)
    {
      gridId0X = extX;
      gridId1X = extX;
      fx = 0;
    }

    if (gridId0Y < 0)
    {
      gridId0Y = 0;
      gridId1Y = 0;
      fy = 0;
    }
    else if (gridId1Y > extY)
    {
      gridId0Y = extY;
      gridId1Y = extY;
      fy = 0;
    }

    if (gridId0Z < 0)
    {
      gridId0Z = 0;
      gridId1Z = 0;
      fz = 0;
    }
    else if (gridId1Z > extZ)
    {
      gridId0Z = extZ;
      gridId1Z = extZ;
      fz = 0;
    }
  }

  // do trilinear interpolation
  vtkIdType factX0 = gridId0X*dataIncs[0];
  vtkIdType factY0 = gridId0Y*dataIncs[1];
  vtkIdType factZ0 = gridId0Z*dataIncs[2];

  vtkIdType factX1 = gridId1X*dataIncs[0];
  vtkIdType factY1 = gridId1Y*dataIncs[1];
  vtkIdType factZ1 = gridId1Z*dataIncs[2];

  T *p000 = dataPtr+factX0+factY0+factZ0;
  T *p001 = dataPtr+factX0+factY0+factZ1;
  T *p010 = dataPtr+factX0+factY1+factZ0;
  T *p011 = dataPtr+factX0+factY1+factZ1;
  T *p100 = dataPtr+factX1+factY0+factZ0;
  T *p101 = dataPtr+factX1+factY0+factZ1;
  T *p110 = dataPtr+factX1+factY1+factZ0;
  T *p111 = dataPtr+factX1+factY1+factZ1;

  float rx = 1.0 - fx;
  float ry = 1.0 - fy;
  float rz = 1.0 - fz;

  float ryrz = ry*rz;
  float ryfz = ry*fz;
  float fyrz = fy*rz;
  float fyfz = fy*fz;

  float rxryrz = rx*ryrz;
  float rxryfz = rx*ryfz;
  float rxfyrz = rx*fyrz;
  float rxfyfz = rx*fyfz;
  float fxryrz = fx*ryrz;
  float fxryfz = fx*ryfz;
  float fxfyrz = fx*fyrz;
  float fxfyfz = fx*fyfz;

  return (T)(rxryrz* *p000 + rxryfz* *p001 + rxfyrz* *p010 + rxfyfz* *p011 +
	     fxryrz* *p100 + fxryfz* *p101 + fxfyrz* *p110 + fxfyfz* *p111);
}

//----------------------------------------------------------------------------
inline static void ComputeHint(int x, int y, int z,
                               float *gridPtr, vtkIdType *gridIncs,
                               int *gridWholeExt,
                               float *inInvSpacing, float *hint)
{
  int extX = gridWholeExt[1]-gridWholeExt[0];
  int extY = gridWholeExt[3]-gridWholeExt[2];
  int extZ = gridWholeExt[5]-gridWholeExt[4];

  int id0X = (x-1) - gridWholeExt[0];
  int id0Y = (y-1) - gridWholeExt[2];
  int id0Z = (z-1) - gridWholeExt[4];

  int id1X = id0X + 2;
  int id1Y = id0Y + 2;
  int id1Z = id0Z + 2;

  if ((id0X | (extX - id1X) |
       id0Y | (extY - id1Y) |
       id0Z | (extZ - id1Z)) < 0)
  {
    hint[0] = hint[1] = hint[2] = 0.0;
    return;
  }
  // do tricubic interpolation
  // this may well be overkill for the computation of the
  // 'hint vector' but the results were markedly better,
  // and the cost is minimal.
  vtkIdType factX0 = id0X*gridIncs[0];
  vtkIdType factY0 = id0Y*gridIncs[1];
  vtkIdType factZ0 = id0Z*gridIncs[2];

  vtkIdType factX = x*gridIncs[0];
  vtkIdType factY = y*gridIncs[1];
  vtkIdType factZ = z*gridIncs[2];

  vtkIdType factX1 = id1X*gridIncs[0];
  vtkIdType factY1 = id1Y*gridIncs[1];
  vtkIdType factZ1 = id1Z*gridIncs[2];

  float *p00 = gridPtr+factX0+factY0+factZ0;
  float *p01 = gridPtr+factX0+factY0+factZ;
  float *p02 = gridPtr+factX0+factY0+factZ1;
  float *p03 = gridPtr+factX0+factY+factZ0;
  float *p04 = gridPtr+factX0+factY+factZ;
  float *p05 = gridPtr+factX0+factY+factZ1;
  float *p06 = gridPtr+factX0+factY1+factZ0;
  float *p07 = gridPtr+factX0+factY1+factZ;
  float *p08 = gridPtr+factX0+factY1+factZ1;

  float *p09 = gridPtr+factX+factY0+factZ0;
  float *p10 = gridPtr+factX+factY0+factZ;
  float *p11 = gridPtr+factX+factY0+factZ1;
  float *p12 = gridPtr+factX+factY+factZ0;
  float *p13 = gridPtr+factX+factY+factZ1;
  float *p14 = gridPtr+factX+factY1+factZ0;
  float *p15 = gridPtr+factX+factY1+factZ;
  float *p16 = gridPtr+factX+factY1+factZ1;

  float *p17 = gridPtr+factX1+factY0+factZ0;
  float *p18 = gridPtr+factX1+factY0+factZ;
  float *p19 = gridPtr+factX1+factY0+factZ1;
  float *p20 = gridPtr+factX1+factY+factZ0;
  float *p21 = gridPtr+factX1+factY+factZ;
  float *p22 = gridPtr+factX1+factY+factZ1;
  float *p23 = gridPtr+factX1+factY1+factZ0;
  float *p24 = gridPtr+factX1+factY1+factZ;
  float *p25 = gridPtr+factX1+factY1+factZ1;

  float *middle = gridPtr+factX+factY+factZ;

  *hint++ = (0.03846 * (*p00++ + *p01++ + *p02++ + *p03++ + *p04++ +
			*p05++ + *p06++ + *p07++ + *p08++ + *p09++ +
			*p10++ + *p11++ + *p12++ + *p13++ + *p14++ +
			*p15++ + *p16++ + *p17++ + *p18++ + *p19++ +
			*p20++ + *p21++ + *p22++ + *p23++ + *p24++ +
			*p25++) - *middle++) * inInvSpacing[0];
  *hint++ = (0.03846 * (*p00++ + *p01++ + *p02++ + *p03++ + *p04++ +
			*p05++ + *p06++ + *p07++ + *p08++ + *p09++ +
			*p10++ + *p11++ + *p12++ + *p13++ + *p14++ +
			*p15++ + *p16++ + *p17++ + *p18++ + *p19++ +
			*p20++ + *p21++ + *p22++ + *p23++ + *p24++ +
			*p25++) - *middle++) * inInvSpacing[1];
  *hint++ = (0.03846 * (*p00++ + *p01++ + *p02++ + *p03++ + *p04++ +
			*p05++ + *p06++ + *p07++ + *p08++ + *p09++ +
			*p10++ + *p11++ + *p12++ + *p13++ + *p14++ +
			*p15++ + *p16++ + *p17++ + *p18++ + *p19++ +
			*p20++ + *p21++ + *p22++ + *p23++ + *p24++ +
			*p25++) - *middle++) * inInvSpacing[2];
}

//----------------------------------------------------------------------------
template <class F, class T>
inline static void CorrelationWorkFunction(T *patientPtr,
					   const vtkIdType *patIncs,
					   const int *patWholeInExt,
					   const int *patext,
					   const unsigned char *modelPoints,
					   const int kernelDiameter,
					   const F *disp,
					   const double sqrt_mod_sum_squared,
					   float &correlation)
{
  unsigned long int aSum=0, topSum=0;

  int idxX, idxY, idxZ;

  idxX = idxY = idxZ = kernelDiameter;


  // we compute the interpolation coefficients and needed datapoints for the
  // first point, and then just shift them around by integer values for the
  // other points in the local extent.
  // change point into integer plus fraction
  F fx,fy,fz;
  int floorX = vtkResliceFloor(disp[0],fx);
  int floorY = vtkResliceFloor(disp[1],fy);
  int floorZ = vtkResliceFloor(disp[2],fz);

  int id0X = floorX - patWholeInExt[0];
  int id0Y = floorY - patWholeInExt[2];
  int id0Z = floorZ - patWholeInExt[4];

  int id1X = id0X + 1;
  int id1Y = id0Y + 1;
  int id1Z = id0Z + 1;

  // do bounds check, most points will be inside so optimize for that
  if ((id0X | (patext[0] - (id1X+kernelDiameter-1)) |
       id0Y | (patext[1] - (id1Y+kernelDiameter-1)) |
       id0Z | (patext[2] - (id1Z+kernelDiameter-1))) < 0)
  {
    correlation = 0.0f;
    return;
  }

  // do trilinear interpolation
  vtkIdType factY0 = id0Y*patIncs[1];
  vtkIdType factZ0 = id0Z*patIncs[2];

  vtkIdType factY1 = id1Y*patIncs[1];
  vtkIdType factZ1 = id1Z*patIncs[2];

  T *p000 = patientPtr+id0X+factY0+factZ0;
  T *p001 = patientPtr+id0X+factY0+factZ1;
  T *p010 = patientPtr+id0X+factY1+factZ0;
  T *p011 = patientPtr+id0X+factY1+factZ1;
  T *p100 = p000 + 1;
  T *p101 = p001 + 1;
  T *p110 = p010 + 1;
  T *p111 = p011 + 1;

  F rx = vtkResliceOneMinusX(fx);
  F ry = vtkResliceOneMinusX(fy);
  F rz = vtkResliceOneMinusX(fz);

  F ryrz = vtkResliceQuikMul(ry*rz);
  F ryfz = vtkResliceQuikMul(ry*fz);
  F fyrz = vtkResliceQuikMul(fy*rz);
  F fyfz = vtkResliceQuikMul(fy*fz);

  F rxryrz = vtkResliceQuikMul(rx*ryrz);
  F rxryfz = vtkResliceQuikMul(rx*ryfz);
  F rxfyrz = vtkResliceQuikMul(rx*fyrz);
  F rxfyfz = vtkResliceQuikMul(rx*fyfz);
  F fxryrz = vtkResliceQuikMul(fx*ryrz);
  F fxryfz = vtkResliceQuikMul(fx*ryfz);
  F fxfyrz = vtkResliceQuikMul(fx*fyrz);
  F fxfyfz = vtkResliceQuikMul(fx*fyfz);

  vtkIdType patContIncY = patIncs[1] - (kernelDiameter);
  vtkIdType patContIncZ = patIncs[2] - (kernelDiameter)*patIncs[1];

  unsigned short patData;

  if (kernelDiameter == 4)  // unroll the loop for a little kernel
  {
    patContIncY = patIncs[1];
    T *p200 = p000 + 2; T *p201 = p001 + 2; T *p210 = p010 + 2; T *p211 = p011 + 2;
    T *p300 = p000 + 3; T *p301 = p001 + 3; T *p310 = p010 + 3; T *p311 = p011 + 3;
    T *p400 = p000 + 4; T *p401 = p001 + 4; T *p410 = p010 + 4; T *p411 = p011 + 4;

    do
    {
      do
      {
	vtkResliceRound(rxryrz* *p000 + rxryfz* *p001 + rxfyrz* *p010 + rxfyfz* *p011 +
  			fxryrz* *p100 + fxryfz* *p101 + fxfyrz* *p110 + fxfyfz* *p111,
  			patData);
	topSum += (unsigned long int) (patData * *modelPoints++);
	aSum += (unsigned long int)(patData * patData);
  	vtkResliceRound(rxryrz* *p100 + rxryfz* *p101 + rxfyrz* *p110 + rxfyfz* *p111 +
  			fxryrz* *p200 + fxryfz* *p201 + fxfyrz* *p210 + fxfyfz* *p211,
  			patData);
	topSum += (unsigned long int) (patData * *modelPoints++);
	aSum += (unsigned long int)(patData * patData);
  	vtkResliceRound(rxryrz* *p200 + rxryfz* *p201 + rxfyrz* *p210 + rxfyfz* *p211 +
  			fxryrz* *p300 + fxryfz* *p301 + fxfyrz* *p310 + fxfyfz* *p311,
  			patData);
	topSum += (unsigned long int) (patData * *modelPoints++);
	aSum += (unsigned long int)(patData * patData);
  	vtkResliceRound(rxryrz* *p300 + rxryfz* *p301 + rxfyrz* *p310 + rxfyfz* *p311 +
  			fxryrz* *p400 + fxryfz* *p401 + fxfyrz* *p410 + fxfyfz* *p411,
  			patData);
	topSum += (unsigned long int) (patData * *modelPoints++);
	aSum += (unsigned long int)(patData * patData);

	p000 += patContIncY; p001 += patContIncY; p010 += patContIncY; p011 += patContIncY;
	p100 += patContIncY; p101 += patContIncY; p110 += patContIncY; p111 += patContIncY;
	p200 += patContIncY; p201 += patContIncY; p210 += patContIncY; p211 += patContIncY;
	p300 += patContIncY; p301 += patContIncY; p310 += patContIncY; p311 += patContIncY;
	p400 += patContIncY; p401 += patContIncY; p410 += patContIncY; p411 += patContIncY;
      }
      while (--idxY);
      p000 += patContIncZ; p001 += patContIncZ; p010 += patContIncZ; p011 += patContIncZ;
      p100 += patContIncZ; p101 += patContIncZ; p110 += patContIncZ; p111 += patContIncZ;
      p200 += patContIncZ; p201 += patContIncZ; p210 += patContIncZ; p211 += patContIncZ;
      p300 += patContIncZ; p301 += patContIncZ; p310 += patContIncZ; p311 += patContIncZ;
      p400 += patContIncZ; p401 += patContIncZ; p410 += patContIncZ; p411 += patContIncZ;
      idxY = kernelDiameter;
    }
    while (--idxZ);
  }

  else // no point unrolling for the other kernels
  {
    do
    {
      do
      {
	do
 {
	  vtkResliceRound(rxryrz* *p000++ + rxryfz* *p001++ +
			  rxfyrz* *p010++ + rxfyfz* *p011++ +
			  fxryrz* *p100++ + fxryfz* *p101++ +
			  fxfyrz* *p110++ + fxfyfz* *p111++,
			  patData);
	  topSum += (unsigned long int) (patData * *modelPoints++);
	  aSum += (unsigned long int)(patData * patData);
 }
	while (--idxX);
	p000 += patContIncY; p001 += patContIncY;
	p010 += patContIncY; p011 += patContIncY;
	p100 += patContIncY; p101 += patContIncY;
	p110 += patContIncY; p111 += patContIncY;
	idxX = kernelDiameter;
      }
      while (--idxY);
      p000 += patContIncZ; p001 += patContIncZ;
      p010 += patContIncZ; p011 += patContIncZ;
      p100 += patContIncZ; p101 += patContIncZ;
      p110 += patContIncZ; p111 += patContIncZ;
      idxY = kernelDiameter;
    }
    while (--idxZ);
  }
  if (sqrt_mod_sum_squared < 0.001 && aSum < 0.00001)
  {
    correlation = 1.0f;
  }
  else
    if (sqrt_mod_sum_squared < 0.001 || aSum < 0.00001)
    {
      correlation = 0.0f;
    }
    else
    {
      correlation = (float)topSum / (sqrt((float)aSum) *
				     sqrt_mod_sum_squared);
    }
}

//----------------------------------------------------------------------------
// stub function since one can not use a templated function in the
// vtkFunctionMinimizer callback
static void CorrelationForExtentAndDisplacement(void *amoebaParmBlock)
{
  _vtkAmoebaParms *pb = (_vtkAmoebaParms *) amoebaParmBlock;

  float disp[3]; // the displacement vector as from amoeba plus the hint
  disp[0] = pb->disps[0] + pb->hint[0];
  disp[1] = pb->disps[1] + pb->hint[1];
  disp[2] = pb->disps[2] + pb->hint[2];
#ifdef VTK_RESLICE_INTEGER_MATH
  int modelDisp[3];
  modelDisp[0] = vtkCastFloatToFixed(disp[0]+ pb->modelExtent[0]);
  modelDisp[1] = vtkCastFloatToFixed(disp[1]+ pb->modelExtent[2]);
  modelDisp[2] = vtkCastFloatToFixed(disp[2]+ pb->modelExtent[4]);
#else
  float modelDisp[3]; // now also add the extent of the model sub-cube
  modelDisp[0] = disp[0] + pb->modelExtent[0];
  modelDisp[1] = disp[1] + pb->modelExtent[2];
  modelDisp[2] = disp[2] + pb->modelExtent[4];
#endif /* VTK_RESLICE_INTEGER_MATH */
  float correlation;

  switch (pb->scalarType)
  {
    vtkTemplateMacro(
      CorrelationWorkFunction((VTK_TT *)(pb->patPtr),
			      pb->patIncs,
			      pb->patWholeInExt,
			      pb->patext,
			      pb->modelPoints,
			      pb->kernelDiameter,
			      modelDisp,
			      pb->sqrt_mod_sum_squared,
			      correlation));
    default:
      cout << "CorrelationForExtentAndDisplacement: Unknown ScalarType\n";
      return;
  }

  float v2 = ((disp[0]*pb->spacing[0])*(disp[0]*pb->spacing[0])+
	      (disp[1]*pb->spacing[1])*(disp[1]*pb->spacing[1])+
	      (disp[2]*pb->spacing[2])*(disp[2]*pb->spacing[2]));
  float v = sqrt(v2);
  float cost;
  v = v*v2;
  if (v < pb->maxLength)
  {
    cost = 0.2 * v / (pb->maxLength - v);
  }
  else
  {
    cost = 1e+38;
  }
  pb->Minimizer->SetScalarResult(cost + 1.0-correlation);
}

//----------------------------------------------------------------------------
// inData[0] is the "patient" data (patData)
// inData[1] is the "model" data   (modData)
// inData[2] is the input grid     (inGrid)
// inData[3] is the stencil, if any
template <class T>
static void vtkImageAmoebaGridExecute(vtkImageAmoebaGrid *self,
				      vtkImageData *patData, T *patPtr,
				      vtkImageData *modData, T *modPtr,
				      vtkImageData *inGrid, float *inGridPtr,
				      vtkImageData *outGrid, float *outPtr,
				      int outExt[6], int &numMinimized,
				      double &vectorLength,
				      double &totalCost, int id)
{
  int idX, idY, idZ;
  int r1, r2, cr1, cr2, iter, rval;
  vtkIdType outIncX, outIncY, outIncZ;

  float *ShrinkFactors = self->GetShrinkFactors();
  int *patWholeInExt = patData->GetWholeExtent();
  int *modWholeInExt = modData->GetWholeExtent();
  int *gridWholeInExt = inGrid->GetWholeExtent();
  vtkIdType *inGridIncs = inGrid->GetIncrements();

  float VectorBounds[3];
  self->GetVectorBounds(VectorBounds);

  vtkFloatingPointType inSpacing[3];
  float inInvSpacing[3];

  float hint[3]; // the mean deformation of the surrounding grid points

  // get the clipping extents
  vtkImageStencilData *stencil = self->GetStencil();

  // since the deformation vectors are calculated in image coordinates, but
  // are applied in world coordinates, they need to be scaled
  // multiply by inInvSpacing to go from world to grid data
  // multiply by spacing to go from grid data to world.
  // Spacing for patData and modData are identical.
  patData->GetSpacing(inSpacing);
  inInvSpacing[0] = 1.0f / inSpacing[0];
  inInvSpacing[1] = 1.0f / inSpacing[1];
  inInvSpacing[2] = 1.0f / inSpacing[2];

  // Get increments to march through output data
  outGrid->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);

  // Set up the parameter block used for the vtkFunctionMinimizer
  _vtkAmoebaParms amoebaParms;
  amoebaParms.patWholeInExt = &patWholeInExt[0];
  amoebaParms.patext[0] = patWholeInExt[1] - patWholeInExt[0];
  amoebaParms.patext[1] = patWholeInExt[3] - patWholeInExt[2];
  amoebaParms.patext[2] = patWholeInExt[5] - patWholeInExt[4];
  amoebaParms.modelData = modData;
  self->GetKernelRadius(amoebaParms.kernelRadius);
  amoebaParms.kernelDiameter = (int)((amoebaParms.kernelRadius[0]+0.5)*2.0);
  int numKernelPix = (amoebaParms.kernelDiameter *
		      amoebaParms.kernelDiameter *
		      amoebaParms.kernelDiameter);
  amoebaParms.modelPoints = new unsigned char[numKernelPix];
  patData->GetIncrements(amoebaParms.patIncs);
  modData->GetIncrements(amoebaParms.modIncs);
  amoebaParms.modelPtr = (void *)modPtr;
  amoebaParms.patPtr = (void *)patPtr;
  amoebaParms.scalarType = patData->GetScalarType();
  amoebaParms.Minimizer = vtkFunctionMinimizer::New();
  amoebaParms.Minimizer->SetFunction(CorrelationForExtentAndDisplacement,
				     &amoebaParms);
  amoebaParms.Minimizer->SetFunctionArgDelete(NULL);
  amoebaParms.Minimizer->SetTolerance(self->GetTolerance());
  amoebaParms.maxLength = 64 * vtkResliceFloor(VectorBounds[0]*
					       VectorBounds[0]*
					       VectorBounds[0]);//(mm)
  amoebaParms.spacing = &inSpacing[0];
  amoebaParms.hint = hint;

  // every displacement vector starts out with (0,0,0) which is the center
  // of the range passed in SetScalarVariableBracket.
  amoebaParms.Minimizer->SetScalarVariableBracket("xdisp",
						  -VectorBounds[0]/2.0,
						  VectorBounds[0]/2.0);
  amoebaParms.Minimizer->SetScalarVariableBracket("ydisp",
						  -VectorBounds[1]/2.0,
						  VectorBounds[1]/2.0);
  amoebaParms.Minimizer->SetScalarVariableBracket("zdisp",
						  -VectorBounds[2]/2.0,
						  VectorBounds[2]/2.0);
  amoebaParms.disps = amoebaParms.Minimizer->GetScalarVarPtr();
  long int iterations = 0;
  numMinimized = 0;
  vectorLength = 0.0f;
  totalCost = 0.0f;
  T temp;

  // Loop through output dataset
  for (idZ = outExt[4]; idZ <= outExt[5]; idZ++)
  {
    amoebaParms.modelVector[2] = idZ*ShrinkFactors[2];
    fprintf(stderr,"%d,",id); // lets us know what thread does what
    for (idY = outExt[2]; !self->AbortExecute && idY <= outExt[3]; idY++)
    {
      amoebaParms.modelVector[1] = idY*ShrinkFactors[1];
      iter = 0;
      cr1 = outExt[0];
      for (;;)
      {
	rval = 0;
	r1 = outExt[1] + 1;
	r2 = outExt[1];
	if (stencil)
 {
	  rval = stencil->GetNextExtent(r1, r2, outExt[0], outExt[1],
					idY, idZ, iter);
 }
	
	cr2 = r1 - 1;
	if (!self->GetReverseStencil())
 {
	  // do unchanged portion
	  for (idX = cr1; idX <= cr2; idX++)
   {
	    *outPtr++ = 0.0;
	    *outPtr++ = 0.0;
	    *outPtr++ = 0.0;
   }
 }
	else
 {
          // do stencil portion
          for (idX = cr1; idX <= cr2; idX++)
          {
	    amoebaParms.modelVector[0] = idX*ShrinkFactors[0];
	    // only compute a vector if we are at > 10% max intensity
	    temp = GetDataAtPoint(modPtr, amoebaParms.modIncs,
				  modWholeInExt,  amoebaParms.modelVector);
	    // since input data scaled (0-255) check for >10% of range
	    // ModelBlock returns 0 if too near edge of volume
	    if ((temp >=1) && (GetModelBlock(modPtr, &amoebaParms)))
     {
	      ComputeHint(idX,idY,idZ,inGridPtr,inGridIncs,
			  gridWholeInExt, inInvSpacing, hint);
              amoebaParms.Minimizer->Minimize();
	
	      // Since the minimization is done in data coordinates,
	      // switch back to world coords which is what the grid is.
	
	      outPtr[0] = inSpacing[0]*(hint[0] + amoebaParms.disps[0]);
	      outPtr[1] = inSpacing[1]*(hint[1] + amoebaParms.disps[1]);
	      outPtr[2] = inSpacing[2]*(hint[2] + amoebaParms.disps[2]);

              vectorLength += sqrt((outPtr[0] * outPtr[0]) +
                                   (outPtr[1] * outPtr[1]) +
                                   (outPtr[2] * outPtr[2]));
	      if (amoebaParms.Minimizer->GetScalarResult() < 1e15)
                totalCost += amoebaParms.Minimizer->GetScalarResult();

              iterations += amoebaParms.Minimizer->GetIterations();
              numMinimized++;
	
	      outPtr+=3;
     }
	    else // target data < 10% of max, no additional displacement
     {
	      *outPtr++ = 0.0;
	      *outPtr++ = 0.0;
	      *outPtr++ = 0.0;
     }
          }
 }
	cr1 = r2 + 1; // for the next kick at the cat
	
	// break if no foreground extents left
	if (rval == 0)
 {
	  break;
 }
	
	if (self->GetReverseStencil())
 {
          // do unchanged portion
          for (idX = r1; idX <= r2; idX++)
          {
	    *outPtr++ = 0.0;
	    *outPtr++ = 0.0;
	    *outPtr++ = 0.0;
          }
 }
	else
 {
	  // do stencil portion
	  for (idX = r1; idX <= r2; idX++)
   {
	    amoebaParms.modelVector[0] = idX*ShrinkFactors[0];
	    // only compute a vector if we are at > 10% max intensity
	    temp = GetDataAtPoint(modPtr, amoebaParms.modIncs,
				  modWholeInExt,  amoebaParms.modelVector);
	    // since input data scaled (0-255) check for >10% of range
	    // ModelBlock returns 0 if too near edge of volume
	    if ((temp >=1) && (GetModelBlock(modPtr, &amoebaParms)))
     {
	      ComputeHint(idX,idY,idZ,inGridPtr,inGridIncs,
			  gridWholeInExt, inInvSpacing, hint);
	      amoebaParms.Minimizer->Minimize();
	
	      // Since the minimization is done in data coordinates,
	      // switch back to world coords which is what the grid is.
		
	      outPtr[0] = inSpacing[0]*(hint[0] + amoebaParms.disps[0]);
	      outPtr[1] = inSpacing[1]*(hint[1] + amoebaParms.disps[1]);
	      outPtr[2] = inSpacing[2]*(hint[2] + amoebaParms.disps[2]);

	      vectorLength += sqrt((outPtr[0] * outPtr[0]) +
              			   (outPtr[1] * outPtr[1]) +
              			   (outPtr[2] * outPtr[2]));

	      if (amoebaParms.Minimizer->GetScalarResult() < 1e15)
                totalCost += amoebaParms.Minimizer->GetScalarResult();

              iterations += amoebaParms.Minimizer->GetIterations();
              numMinimized++;
		
	      outPtr+=3;
     }
	    else // target data < 10% of max, no additional displacement
     {
	      *outPtr++ = 0.0;
	      *outPtr++ = 0.0;
	      *outPtr++ = 0.0;
     }
   }
 }
      }
      outPtr += outIncY;
    }
    outPtr += outIncZ;
  }
  amoebaParms.Minimizer->Delete();
  delete [] amoebaParms.modelPoints;

  if (numMinimized)
  {
    cout << "\nMinimized "<< numMinimized <<" using average "
	 <<iterations/numMinimized<<" funks in thread "<<id<<".\n";
  }
  else
  {
    cout << "\nMinimized 0 in thread "<<id<<".\n";
  }

}

//----------------------------------------------------------------------------
// This method is passed a input and output region, and executes the filter
// algorithm to fill the output from the input.
void vtkImageAmoebaGrid::ThreadedExecute(vtkImageData **inData,
					 vtkImageData *outData,
					 int outExt[6], int id)
{
  void  *patPtr, *modPtr, *inGridPtr, *outPtr;

  vtkDebugMacro(<< "Execute: inData = " << inData
                << ", outData = " << outData);

  if (inData[0] == NULL)
  {
    vtkErrorMacro(<< "Input 0 must be specified (patient data).");
    return;
  }
  patPtr = inData[0]->GetScalarPointer();

  if (inData[1] == NULL)
  {
    vtkErrorMacro(<< "Input 1 must be specified (model data).");
    return;
  }
  modPtr = inData[1]->GetScalarPointer();

  if (inData[2] == NULL)
  {
    vtkErrorMacro(<< "Input 2 must be specified (input grid).");
    return;
  }
  inGridPtr = inData[2]->GetScalarPointer();

  outPtr = outData->GetScalarPointerForExtent(outExt);

  // this filter expects that the image inputs have the same number of components
  if (inData[0]->GetNumberOfScalarComponents() !=
      inData[1]->GetNumberOfScalarComponents())
  {
    vtkErrorMacro(<< "Execute: input1 NumberOfScalarComponents, "
                  << inData[0]->GetNumberOfScalarComponents()
                  << ", must match out input2 NumberOfScalarComponents "
                  << inData[1]->GetNumberOfScalarComponents());
    return;
  }

  // This filter expects that input grid is the same type as output,
  // which is set to VTK_FLOAT.
  if (inData[2]->GetScalarType() != outData->GetScalarType())
  {
    vtkErrorMacro("Execute: Grid Input ScalarType, " << inData[2]->GetScalarType()
		  << ", must match Output ScalarType (VTK_FLOAT)"
		  << outData->GetScalarType());
    return;
  }

  switch (inData[0]->GetScalarType())
  {
    case VTK_CHAR:
      vtkImageAmoebaGridExecute(this,
				inData[0], (char *)patPtr,
				inData[1], (char *)modPtr,
				inData[2], (float *)inGridPtr,
				outData,(float *)(outPtr),
				outExt, this->VectorsMinimized[id],
				this->VectorLength[id],
				this->TotalCost[id], id);
      break;
    case VTK_UNSIGNED_CHAR:
      vtkImageAmoebaGridExecute(this,
				inData[0], (unsigned char *)patPtr,
				inData[1], (unsigned char *)modPtr,
				inData[2], (float *)inGridPtr,
				outData,(float *)(outPtr),
				outExt, this->VectorsMinimized[id],
				this->VectorLength[id],
				this->TotalCost[id], id);
      break;
    case VTK_SHORT:
      vtkImageAmoebaGridExecute(this,
				inData[0], (short *)patPtr,
				inData[1], (short *)modPtr,
				inData[2], (float *)inGridPtr,
				outData,(float *)(outPtr),
				outExt, this->VectorsMinimized[id],
				this->VectorLength[id],
				this->TotalCost[id], id);
      break;
    case VTK_UNSIGNED_SHORT:
      vtkImageAmoebaGridExecute(this,
				inData[0], (unsigned short *)patPtr,
				inData[1], (unsigned short *)modPtr,
				inData[2], (float *)inGridPtr,
				outData,(float *)(outPtr),
				outExt, this->VectorsMinimized[id],
				this->VectorLength[id],
				this->TotalCost[id], id);
      break;
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
      return;
  }
}

//----------------------------------------------------------------------------
void vtkImageAmoebaGrid::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkImageMultipleInputFilter::PrintSelf(os,indent);

  os << indent << "Stencil: " << this->GetStencil() << "\n";
  os << indent << "ReverseStencil: " << (this->ReverseStencil ?
		                         "On\n" : "Off\n");

  os << indent << "ShrinkFactors: (" << this->ShrinkFactors[0] << ", "
     << this->ShrinkFactors[1] << ", " << this->ShrinkFactors[2] << ")\n";
  os << indent << "KernelRadius: (" << this->KernelRadius[0] << ", "
     << this->KernelRadius[1] << ", " << this->KernelRadius[2] << ")\n";

}

