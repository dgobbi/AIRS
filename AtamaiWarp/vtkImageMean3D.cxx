/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImageMean3D.cxx,v $
  Language:  C++
  Date:      $Date: 2007/08/24 20:02:25 $
  Version:   $Revision: 1.7 $
  Thanks:    Thanks to C. Charles Law who developed this class.

Copyright (c) 1993-2001 Ken Martin, Will Schroeder, Bill Lorensen 
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
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#include "vtkImageMean3D.h"
#include "vtkImageData.h"
#include "vtkObjectFactory.h"

//----------------------------------------------------------------------------
vtkImageMean3D* vtkImageMean3D::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkImageMean3D");
  if(ret)
    {
    return (vtkImageMean3D*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkImageMean3D;
}

//----------------------------------------------------------------------------
// Construct an instance of vtkImageMean3D fitler.
vtkImageMean3D::vtkImageMean3D()
{
  this->SetKernelSize(1,1,1);
  this->HandleBoundaries = 1;
  this->SmoothThreshold = 0.0;
}

//----------------------------------------------------------------------------
void vtkImageMean3D::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "NumberOfElements: " << this->NumberOfElements << endl;
  os << indent << "SmoothThreshold: " << this->SmoothThreshold << endl;
}

//----------------------------------------------------------------------------
// This method sets the size of the neighborhood.  It also sets the 
// default middle of the neighborhood 
void vtkImageMean3D::SetKernelSize(int size0, int size1, int size2)
{  
  int volume;
  int modified = 1;
  
  if (this->KernelSize[0] == size0 && this->KernelSize[1] == size1 && 
      this->KernelSize[2] == size2)
    {
    modified = 0;
    }
  
  // Set the kernel size and middle
  volume = 1;
  this->KernelSize[0] = size0;
  this->KernelMiddle[0] = size0 / 2;
  volume *= size0;
  this->KernelSize[1] = size1;
  this->KernelMiddle[1] = size1 / 2;
  volume *= size1;
  this->KernelSize[2] = size2;
  this->KernelMiddle[2] = size2 / 2;
  volume *= size2;

  this->NumberOfElements = volume;
  if ( modified )
  {
  this->Modified();
  }
}

//----------------------------------------------------------------------------
// This method contains the second switch statement that calls the correct
// templated function for the mask types.
template <class T>
static void vtkImageMean3DExecute(vtkImageMean3D *self,
				  vtkImageData *inData, T *inPtr, 
				  vtkImageData *outData, T *outPtr,
				  int outExt[6], int id)
{
  int *kernelMiddle, *kernelSize;
  int NumberOfElements;
  unsigned long count = 0; // keep count for the progress meter
  unsigned long target;
  // For looping though output (and input) pixels.
  int outIdx0, outIdx1, outIdx2;
  vtkIdType inInc0, inInc1, inInc2;
  vtkIdType outIncX, outIncY, outIncZ;
  T *inPtr0, *inPtr1, *inPtr2;
  // For looping through hood pixels
  int hoodMin0, hoodMax0, hoodMin1, hoodMax1, hoodMin2, hoodMax2;
  int hoodStartMin0, hoodStartMax0, hoodStartMin1, hoodStartMax1;
  int hoodIdx0, hoodIdx1, hoodIdx2;
  T *tmpPtr0, *tmpPtr1, *tmpPtr2;

  // The portion of the out image that needs no boundary processing.
  int middleMin0, middleMax0, middleMin1, middleMax1, middleMin2, middleMax2;
  int numComp;
  // variables for the mean calc
  int *inExt;
  double sumX, sumY, sumZ;
  int numPixels;
  T centerPixelX, centerPixelY, centerPixelZ;
  float smoothThreshold = self->GetSmoothThreshold();
  float centerWeighting = self->GetCenterWeighting();
  float surroundWeighting = self->GetSurroundWeighting();

  // avoid compiler warnings
  centerPixelX = centerPixelY = centerPixelZ = (T)0;
  
  // Get information to march through data
  inData->GetIncrements(inInc0, inInc1, inInc2); 
  outData->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);
  kernelMiddle = self->GetKernelMiddle();
  kernelSize = self->GetKernelSize();

  numComp = inData->GetNumberOfScalarComponents();

  hoodMin0 = outExt[0] - kernelMiddle[0]; 
  hoodMin1 = outExt[2] - kernelMiddle[1]; 
  hoodMin2 = outExt[4] - kernelMiddle[2]; 
  hoodMax0 = kernelSize[0] + hoodMin0 - 1;
  hoodMax1 = kernelSize[1] + hoodMin1 - 1;
  hoodMax2 = kernelSize[2] + hoodMin2 - 1;
  
  // Clip by the input image extent
  inExt = inData->GetExtent();
  hoodMin0 = (hoodMin0 > inExt[0]) ? hoodMin0 : inExt[0];
  hoodMin1 = (hoodMin1 > inExt[2]) ? hoodMin1 : inExt[2];
  hoodMin2 = (hoodMin2 > inExt[4]) ? hoodMin2 : inExt[4];
  hoodMax0 = (hoodMax0 < inExt[1]) ? hoodMax0 : inExt[1];
  hoodMax1 = (hoodMax1 < inExt[3]) ? hoodMax1 : inExt[3];
  hoodMax2 = (hoodMax2 < inExt[5]) ? hoodMax2 : inExt[5];

  // Save the starting neighborhood dimensions (2 loops only once)
  hoodStartMin0 = hoodMin0;    hoodStartMax0 = hoodMax0;
  hoodStartMin1 = hoodMin1;    hoodStartMax1 = hoodMax1;
  
  // The portion of the output that needs no boundary computation.
  middleMin0 = inExt[0] + kernelMiddle[0];
  middleMax0 = inExt[1] - (kernelSize[0] - 1) + kernelMiddle[0];
  middleMin1 = inExt[2] + kernelMiddle[1];
  middleMax1 = inExt[3] - (kernelSize[1] - 1) + kernelMiddle[1];
  middleMin2 = inExt[4] + kernelMiddle[2];
  middleMax2 = inExt[5] - (kernelSize[2] - 1) + kernelMiddle[2];
  
  target = (unsigned long)((outExt[5] - outExt[4] + 1)*
			   (outExt[3] - outExt[2] + 1)/50.0);
  target++;
  
  NumberOfElements = self->GetNumberOfElements();
  
  long int numSmoothed=0;
  long int numIgnored = 0;
  
  // loop through pixel of output
  inPtr = (T *)inData->GetScalarPointer(hoodMin0,hoodMin1,hoodMin2);
  inPtr2 = inPtr;
  for (outIdx2 = outExt[4]; outIdx2 <= outExt[5]; ++outIdx2)
    {
    inPtr1 = inPtr2;
    hoodMin1 = hoodStartMin1;
    hoodMax1 = hoodStartMax1;
    for (outIdx1 = outExt[2]; 
	 !self->AbortExecute && outIdx1 <= outExt[3]; ++outIdx1)
      {
      if (!id) 
	{
	if (!(count%target))
	  {
	  self->UpdateProgress(count/(50.0*target));
	  }
	count++;
	}
      inPtr0 = inPtr1;
      hoodMin0 = hoodStartMin0;
      hoodMax0 = hoodStartMax0;
      for (outIdx0 = outExt[0]; outIdx0 <= outExt[1]; ++outIdx0)
	{
	// Compute mean of neighborhood
	// loop through neighborhood pixels
	// This could be sped up pretty dramatically by only adding the 
        // pixels on the leading edge of the neighborhood and subtracting
        // those on the trailing edge of the neighborhood as it moves
        // through the input data.
	tmpPtr2 = inPtr0;
	  
	sumX = sumY = sumZ = 0.0;
	numPixels = 0;
	for (hoodIdx2 = hoodMin2; hoodIdx2 <= hoodMax2; ++hoodIdx2)
	  {
	  tmpPtr1 = tmpPtr2;
	  for (hoodIdx1 = hoodMin1; hoodIdx1 <= hoodMax1; ++hoodIdx1)
	    {
	    tmpPtr0 = tmpPtr1;
	    for (hoodIdx0 = hoodMin0; hoodIdx0 <= hoodMax0; ++hoodIdx0)
	      {
	      // Add this pixel to the mean
	      if ((outIdx0-hoodIdx0) | (outIdx1-hoodIdx1) | (outIdx2-hoodIdx2))
		{
		numPixels++;
		sumX+=(double)*(tmpPtr0++);
		sumY+=(double)*(tmpPtr0++);
		sumZ+=(double)*(tmpPtr0++);
		}
	      else // unless it is the center pixel, in which case we save it.
		{
		centerPixelX = *(tmpPtr0++);
		centerPixelY = *(tmpPtr0++);
		centerPixelZ = *(tmpPtr0++);
		}
	      }
	    tmpPtr1 += inInc1;
	    }
	  tmpPtr2 += inInc2;
	  }
	
	// Replace this vector with the hood mean assuming its longer 
	// than the threshold length.
	if (((centerPixelX*centerPixelX + 
	      centerPixelY*centerPixelY + 
	      centerPixelZ*centerPixelZ) > smoothThreshold) &&(numPixels))
	  {
	  numSmoothed++;
	  *(outPtr++) = (T)(centerWeighting*centerPixelX + 
			    surroundWeighting*(sumX/numPixels));
	  *(outPtr++) = (T)(centerWeighting*centerPixelY + 
			    surroundWeighting*(sumY/numPixels));
	  *(outPtr++) = (T)(centerWeighting*centerPixelZ + 
			    surroundWeighting*(sumZ/numPixels));
	  } 
	else 
	  {
	  numIgnored++;
	  *(outPtr++) = (T)centerPixelX;
	  *(outPtr++) = (T)centerPixelY;
	  *(outPtr++) = (T)centerPixelZ;
	  }	    
	  
	// shift neighborhood considering boundaries
	if (outIdx0 >= middleMin0)
	  {
	  inPtr0 += inInc0;
	  ++hoodMin0;
	  }
	if (outIdx0 < middleMax0)
	  {
	  ++hoodMax0;
	  }
	}
      // shift neighborhood considering boundaries
      if (outIdx1 >= middleMin1)
	{
	inPtr1 += inInc1;
	++hoodMin1;
	}
      if (outIdx1 < middleMax1)
	{
	++hoodMax1;
	}
      outPtr += outIncY;
      }
    // shift neighborhood considering boundaries
    if (outIdx2 >= middleMin2)
      {
      inPtr2 += inInc2;
      ++hoodMin2;
      }
    if (outIdx2 < middleMax2)
      {
      ++hoodMax2;
      }
    outPtr += outIncZ;
    }
  
  cout << "Smoothed "<<numSmoothed<<". Ignored "<<numIgnored<<".\n";
}

//----------------------------------------------------------------------------
// This method contains the first switch statement that calls the correct
// templated function for the input and output region types.
void vtkImageMean3D::ThreadedExecute(vtkImageData *inData, 
				       vtkImageData *outData,
				       int outExt[6], int id)
{
  void *inPtr = inData->GetScalarPointerForExtent(outExt);
  void *outPtr = outData->GetScalarPointerForExtent(outExt);
  
  vtkDebugMacro(<< "Execute: inData = " << inData 
  << ", outData = " << outData);
  
  // this filter expects that input is the same type as output.
  if (inData->GetScalarType() != outData->GetScalarType())
    {
    vtkErrorMacro(<< "Execute: input ScalarType, " << inData->GetScalarType()
      << ", must match out ScalarType " << outData->GetScalarType());
    return;
    }
  
  switch (inData->GetScalarType())
    {
#if (VTK_MAJOR_VERSION < 5)
    vtkTemplateMacro7(vtkImageMean3DExecute, this,inData, (VTK_TT *)(inPtr), 
                      outData, (VTK_TT *)(outPtr),outExt, id);
#else
    vtkTemplateMacro(
        vtkImageMean3DExecute(this,inData, (VTK_TT *)(inPtr), 
                      outData, (VTK_TT *)(outPtr),outExt, id));
#endif
    default:
      vtkErrorMacro(<< "Execute: Unknown input ScalarType");
      return;
    }
}


