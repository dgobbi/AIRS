/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImageCrossCorrelation.cxx,v $
  Language:  C++
  Date:      $Date: 2004/07/13 14:43:11 $
  Version:   $Revision: 1.1 $

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
#include "vtkImageData.h"
#include "vtkImageCrossCorrelation.h"
#include "vtkObjectFactory.h"



//------------------------------------------------------------------------------
vtkImageCrossCorrelation* vtkImageCrossCorrelation::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkImageCrossCorrelation");
  if(ret)
    {
    return (vtkImageCrossCorrelation*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkImageCrossCorrelation;
}


//----------------------------------------------------------------------------
// fast floor() function for converting a float to an int
// (the floor() implementation on some computers is much slower than this,
// because they require some 'exact' behaviour that we don't).
// convert a float into an integer plus a fraction  
static inline int vtkResliceFloor(double x)
{
#if defined mips || defined sparc
  return (int)((unsigned int)(x + 2147483648.0) - 2147483648U);
#elif defined i386
  unsigned int hilo[2];
  *((double *)hilo) = x + 103079215104.0;  // (2**(52-radix))*1.5
  return (int)((hilo[1]<<16)|(hilo[0]>>16));
#else
  return int(floor(x));
#endif
}

static inline int vtkResliceCeil(double x)
{
  return -vtkResliceFloor(-x - 1.0) - 1;
}

static inline int vtkResliceFloor(float x, float &f)
{
  int ix = vtkResliceFloor(x);
  f = x - ix;

  return ix;
}





//----------------------------------------------------------------------------
vtkImageCrossCorrelation::vtkImageCrossCorrelation()
{
  this->ShrinkFactors[0] = this->ShrinkFactors[1] = this->ShrinkFactors[2] = 1.0f;
  this->KernelRadius[0] = this->KernelRadius[1] = this->KernelRadius[2] = 1.0f;
}


//----------------------------------------------------------------------------
// This method computes the Region of input necessary to generate outRegion.
void vtkImageCrossCorrelation::ComputeInputUpdateExtent(int inExt[6], 
							int outExt[6],
							int vtkNotUsed(whichInput))
{
  int idx;

  int *wholeInExt = this->GetInput()->GetWholeExtent();

  for (idx = 0; idx < 3; ++idx)
    {
      
    // For Min.
    inExt[idx*2] = vtkResliceFloor(outExt[idx*2]*this->ShrinkFactors[idx]-
				   this->KernelRadius[idx]);
    if (inExt[idx*2] < wholeInExt[idx*2])
      {
      inExt[idx*2] = wholeInExt[idx*2];
      }

    // For Max.
    inExt[idx*2+1] = vtkResliceCeil(outExt[idx*2+1]*this->ShrinkFactors[idx]+
				    this->KernelRadius[idx]);
    if (inExt[idx*2+1] > wholeInExt[idx*2+1]) 
      { 
      inExt[idx*2+1] = wholeInExt[idx*2+1];
      }
    }
}




//----------------------------------------------------------------------------
//Set Up the output
//Force 1 component floats
// Set spacing so that the bounds of the output mimic those of the input
// Bounds extend to the edges of the volume, but no cross correlation info
// is calculated there. 
void vtkImageCrossCorrelation::ExecuteInformation(vtkImageData **inDatas, 
						  vtkImageData *outData)
{
  int idx;
  int wholeInExt[6];
  double spacing[3];
  
  inDatas[0]->Update();
  inDatas[1]->Update();
  inDatas[0]->GetWholeExtent(wholeInExt);
  inDatas[0]->GetSpacing(spacing);


  for (idx = 0; idx < 3; ++idx)
    {
      // Scale the output extent
      this->ShrinkFactors[idx] = 
	((float)wholeInExt[2*idx+1] - (float)wholeInExt[2*idx])/ 
	((float)this->OutputExtent[2*idx+1] - (float)this->OutputExtent[2*idx]);
      // Change the data spacing
      spacing[idx] *= this->ShrinkFactors[idx];
    }

  outData->SetWholeExtent(this->OutputExtent);
  outData->SetSpacing(spacing);

  outData->SetOrigin(inDatas[0]->GetOrigin());
  outData->SetNumberOfScalarComponents(1);
  outData->SetScalarType(VTK_FLOAT);
  
}




//----------------------------------------------------------------------------
// Ass-you-me-ptions:
// input data match in extent, spacing
// 1 scalar component only
// all 3 kernel dimensions are equal!
template <class T>
static void vtkImageCrossCorrelationExecute(vtkImageCrossCorrelation *self, 
					    T *in1Ptr,
					    T *in2Ptr,
					    vtkImageData **inData,
					    vtkImageData *outData,
					    float *outPtr, int outExt[6],
					    int id)
{
  int outIdxX, outIdxY, outIdxZ;
  int outIncX, outIncY, outIncZ;

  int inIdxX, inIdxY, inIdxZ;
  double aSum = 0.0;
  double bSum = 0.0;
  double topSum = 0.0;

  float *ShrinkFactors;
  float *KernelRadius;
  int *wholeInExt;

  float fx,fy,fz;
  float rx,ry,rz;
  float ryrz, ryfz, fyrz, fyfz;
  int floorX, floorY, floorZ;
  int id0X, id0Y, id0Z;
  int id1X, id1Y, id1Z;

  int irxryrz, irxryfz, irxfyrz, irxfyfz;
  int ifxryrz, ifxryfz, ifxfyrz, ifxfyfz;

  int factY0, factY1, factZ0, factZ1;

  wholeInExt = inData[0]->GetWholeExtent();
  int extX = wholeInExt[1] - wholeInExt[0];
  int extY = wholeInExt[3] - wholeInExt[2];
  int extZ = wholeInExt[5] - wholeInExt[4];

  T *pA000; T *pA001; T *pA010; T *pA011; // interpolation pointers for imgA
  T *pA100; T *pA101; T *pA110; T *pA111;

  T *pB000; T *pB001; T *pB010; T *pB011; // interpolation pointers for imgB
  T *pB100; T *pB101; T *pB110; T *pB111;

  long int data1,data2; // temp holders for the interpolated result

  ShrinkFactors = self->GetShrinkFactors();
  KernelRadius = self->GetKernelRadius();

  int kernelDiameter = (int)((KernelRadius[0]+0.5)*2.0);
  int inIncX, inIncY, inIncZ;
  inData[0]->GetIncrements(inIncX, inIncY, inIncZ);
  int contInIncY = inIncY - (kernelDiameter);
  int contInIncZ = inIncZ - (kernelDiameter)*inIncY;

  // Get increments to march through data 
  outData->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);

  // Loop through output dataset
  for (outIdxZ = outExt[4]; outIdxZ <= outExt[5]; outIdxZ++)
    {
    floorZ = vtkResliceFloor((outIdxZ*ShrinkFactors[2]) - KernelRadius[2],fz);
    id0Z = floorZ - wholeInExt[4];
    id1Z = floorZ+1;
    factZ0 = id0Z * inIncZ;
    factZ1 = id1Z * inIncZ;
    rz = 1.0 - fz;
   
    for (outIdxY = outExt[2]; 
	 !self->AbortExecute && outIdxY <= outExt[3]; outIdxY++)
      {
      floorY = vtkResliceFloor((outIdxY*ShrinkFactors[1]) 
			       - KernelRadius[1],fy);
      id0Y = floorY - wholeInExt[2];
      id1Y = floorY+1;
      factY0 = id0Y * inIncY;
      factY1 = id1Y * inIncY;
      ry = 1.0 - fy;
      ryrz = ry*rz;
      ryfz = ry*fz;
      fyrz = fy*rz;
      fyfz = fy*fz;

      for (outIdxX = outExt[0]; outIdxX <= outExt[1]; outIdxX++)
	{
	floorX = vtkResliceFloor((outIdxX*ShrinkFactors[0]) 
				 - KernelRadius[0],fx);
	id0X = floorX - wholeInExt[0];
	id1X = floorX+1;
	if ((id0X | (extX - (id1X+kernelDiameter-1)) |
	     id0Y | (extY - (id1Y+kernelDiameter-1)) |
	     id0Z | (extZ - (id1Z+kernelDiameter-1))) < 0)
	  {
	  *outPtr++ = 1.0f; // edges are perfect correlation!
	  }
	else
	  {
	  rx = 1.0-fx;
	  pA000=in1Ptr+id0X+factY0+factZ0;    pB000=in2Ptr+id0X+factY0+factZ0;
	  pA001=in1Ptr+id0X+factY0+factZ1;    pB001=in2Ptr+id0X+factY0+factZ1;
	  pA010=in1Ptr+id0X+factY1+factZ0;    pB010=in2Ptr+id0X+factY1+factZ0;
	  pA011=in1Ptr+id0X+factY1+factZ1;    pB011=in2Ptr+id0X+factY1+factZ1;
	  pA100=pA000+1;    pA101=pA001+1;    pA110=pA010+1;    pA111=pA011+1;
	  pB100=pB000+1;    pB101=pB001+1;    pB110=pB010+1;    pB111=pB011+1;

	  irxryrz = (int)(rx*ryrz *512);
	  irxryfz = (int)(rx*ryfz *512);
	  irxfyrz = (int)(rx*fyrz *512);
	  irxfyfz = (int)(rx*fyfz *512);
	  ifxryrz = (int)(fx*ryrz *512);
	  ifxryfz = (int)(fx*ryfz *512);
	  ifxfyrz = (int)(fx*fyrz *512);
	  ifxfyfz = (int)(fx*fyfz *512);

	  inIdxX = kernelDiameter;
	  inIdxY = kernelDiameter;
	  inIdxZ = kernelDiameter;

	  aSum = 0.0;
	  bSum = 0.0;
	  topSum = 0.0;

	  do
	    {
	    do
	      {
	      do
		{
		data1 = ((long int) (irxryrz* *pA000++ + irxryfz* *pA001++ + 
				     irxfyrz* *pA010++ + irxfyfz* *pA011++ +
				     ifxryrz* *pA100++ + ifxryfz* *pA101++ + 
				     ifxfyrz* *pA110++ + ifxfyfz* *pA111++)) >> 9;
		data2 = ((long int) (irxryrz* *pB000++ + irxryfz* *pB001++ + 
				     irxfyrz* *pB010++ + irxfyfz* *pB011++ +
				     ifxryrz* *pB100++ + ifxryfz* *pB101++ + 
				     ifxfyrz* *pB110++ + ifxfyfz* *pB111++)) >> 9;
		
		topSum +=  fabs((double)data1 * (double)data2);
		aSum += ((double)data1 * (double)data1);
		bSum += ((double)data2 * (double)data2);
		}
	      while (--inIdxX);
	      pA000+=contInIncY;pA001+=contInIncY;pA010+=contInIncY;pA011+=contInIncY;
	      pA100+=contInIncY;pA101+=contInIncY;pA110+=contInIncY;pA111+=contInIncY;
	      pB000+=contInIncY;pB001+=contInIncY;pB010+=contInIncY;pB011+=contInIncY;
	      pB100+=contInIncY;pB101+=contInIncY;pB110+=contInIncY;pB111+=contInIncY;
	      inIdxX = kernelDiameter;
	      }
	    while (--inIdxY);
	    pA000+=contInIncZ;pA001+=contInIncZ;pA010+=contInIncZ;pA011+=contInIncZ;
	    pA100+=contInIncZ;pA101+=contInIncZ;pA110+=contInIncZ;pA111+=contInIncZ;
	    pB000+=contInIncZ;pB001+=contInIncZ;pB010+=contInIncZ;pB011+=contInIncZ;
	    pB100+=contInIncZ;pB101+=contInIncZ;pB110+=contInIncZ;pB111+=contInIncZ;
	    inIdxY = kernelDiameter;
	    }
	  while (--inIdxZ);
        
	  if ((aSum>0) && (bSum>0))
	    {
	    *outPtr++ = topSum / (sqrt(aSum) * sqrt(bSum));
	    }
	  else
	    {
	    *outPtr++ = 1.0f; // All zeros counts as perfect correlation
	    }
	  }
	}
      outPtr += outIncY;
      }
    outPtr += outIncZ;
    }
}



//----------------------------------------------------------------------------
// This method is passed a input and output region, and executes the filter
// algorithm to fill the output from the input.
// It just executes a switch statement to call the correct function for
// the regions data types.
void vtkImageCrossCorrelation::ThreadedExecute(vtkImageData **inData, 
					       vtkImageData *outData,
					       int outExt[6], int id) 
{
  void *inPtr1;
  void *outPtr;
  int wholeInExt[6];
  int wholeInExt2[6];

  vtkDebugMacro(<< "Execute: inData = " << inData 
		<< ", outData = " << outData);

  if (inData[0] == NULL)
    {
    vtkErrorMacro(<< "Input 0 must be specified.");
    return;
    }

  inData[0]->GetWholeExtent(wholeInExt);

  inPtr1 = inData[0]->GetScalarPointer();
  
  void *inPtr2;

  if (inData[1] == NULL)
    {
    vtkErrorMacro(<< "Input " << 1 << " must be specified.");
    return;
    }

  inData[1]->GetWholeExtent(wholeInExt2);
  if ((wholeInExt[0] - wholeInExt2[0]) | (wholeInExt[1] - wholeInExt2[1]) |
      (wholeInExt[2] - wholeInExt2[2]) | (wholeInExt[3] - wholeInExt2[3]) |
      (wholeInExt[4] - wholeInExt2[4]) | (wholeInExt[5] - wholeInExt2[5]))
    {
    vtkErrorMacro(<<"Inputs don't have matching wholeExtents.");
    return;
    }
  inPtr2 = inData[1]->GetScalarPointer();
  outPtr = outData->GetScalarPointerForExtent(outExt);


  
  // this filter expects that inputs that have the same number of components
  if (inData[0]->GetNumberOfScalarComponents() != 
      inData[1]->GetNumberOfScalarComponents())
    {
    vtkErrorMacro(<< "Execute: input1 NumberOfScalarComponents, "
    << inData[0]->GetNumberOfScalarComponents()
    << ", must match out input2 NumberOfScalarComponents "
    << inData[1]->GetNumberOfScalarComponents());
    return;
    }

  // scalar types of the input must be the same
  if (inData[0]->GetScalarType() != inData[1]->GetScalarType())
    {
    vtkErrorMacro(<< "Execute: ScalarType Mismatch");
    return;
    }
  
  switch (inData[0]->GetScalarType())
    {
    vtkTemplateMacro8(vtkImageCrossCorrelationExecute,
		      this, (VTK_TT *)(inPtr1), (VTK_TT *)(inPtr2),inData,outData,
		      (float *)(outPtr), outExt,
		      id);
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
      return;
    }
}


void vtkImageCrossCorrelation::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkImageTwoInputFilter::PrintSelf(os,indent);

  os << indent << "ShrinkFactors: (" << this->ShrinkFactors[0] << ", "
     << this->ShrinkFactors[1] << ", " << this->ShrinkFactors[2] << ")\n";
  os << indent << "KernelRadius: (" << this->KernelRadius[0] << ", "
     << this->KernelRadius[1] << ", " << this->KernelRadius[2] << ")\n";

}

