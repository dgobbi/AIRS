/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImageMutualInformation.cxx,v $
  Language:  C++
  Date:      $Date: 2006/09/21 13:30:37 $
  Version:   $Revision: 1.7 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageMutualInformation.h"

#include "vtkImageData.h"
#include "vtkImageStencilData.h"
#include "vtkObjectFactory.h"

#include <math.h>

vtkCxxRevisionMacro(vtkImageMutualInformation, "$Revision: 1.7 $");
vtkStandardNewMacro(vtkImageMutualInformation);

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageMutualInformation::vtkImageMutualInformation()
{
  this->ImageAComponentSpacing = 1.0;
  this->ImageAComponentOrigin = 0.0;
  this->ImageAComponentExtent[0] = 0;
  this->ImageAComponentExtent[1] = 255;
  this->ImageBComponentSpacing = 1.0;
  this->ImageBComponentOrigin = 0.0;
  this->ImageBComponentExtent[0] = 0;
  this->ImageBComponentExtent[1] = 255;
  
  this->ReverseStencil = 0;

  this->NormalizedMI = 0.0;
}


//----------------------------------------------------------------------------
vtkImageMutualInformation::~vtkImageMutualInformation()
{
}

//----------------------------------------------------------------------------
void vtkImageMutualInformation::SetImageAComponentExtent(int extent[2])
{
  if (this->ImageAComponentExtent[0] != extent[0])
    {
    this->ImageAComponentExtent[0] = extent[0];
    this->Modified();
    }
  if (this->ImageAComponentExtent[1] != extent[1])
    {
    this->ImageAComponentExtent[1] = extent[1];
    this->Modified();
    }
}
//----------------------------------------------------------------------------
void vtkImageMutualInformation::SetImageBComponentExtent(int extent[2])
{
  if (this->ImageBComponentExtent[0] != extent[0])
    {
    this->ImageBComponentExtent[0] = extent[0];
    this->Modified();
    }
  if (this->ImageBComponentExtent[1] != extent[1])
    {
    this->ImageBComponentExtent[1] = extent[1];
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkImageMutualInformation::SetImageAComponentExtent(int min, int max)
{
  int extent[2];
  
  extent[0] = min;  extent[1] = max;
  this->SetImageAComponentExtent(extent);
}
//----------------------------------------------------------------------------
void vtkImageMutualInformation::SetImageBComponentExtent(int min, int max)
{
  int extent[2];
  
  extent[0] = min;  extent[1] = max;
  this->SetImageBComponentExtent(extent);
}


//----------------------------------------------------------------------------
void vtkImageMutualInformation::GetImageAComponentExtent(int extent[2])
{
  extent[0] = this->ImageAComponentExtent[0];
  extent[1] = this->ImageAComponentExtent[1];
}
//----------------------------------------------------------------------------
void vtkImageMutualInformation::GetImageBComponentExtent(int extent[2])
{
  extent[0] = this->ImageBComponentExtent[0];
  extent[1] = this->ImageBComponentExtent[1];
}


//----------------------------------------------------------------------------
void vtkImageMutualInformation::SetStencil(vtkImageStencilData *stencil)
{
  this->vtkProcessObject::SetNthInput(2, stencil); 
}


//----------------------------------------------------------------------------
vtkImageStencilData *vtkImageMutualInformation::GetStencil()
{
  if (this->NumberOfInputs < 3) 
    { 
    return NULL;
    }
  else
    {
    return (vtkImageStencilData *)(this->Inputs[2]); 
    }
}


//--------------------------------------------------------------------------
// The 'floor' function on x86, mips and ppc is many times slower than these
// and is used a lot in this code, optimize for different CPU architectures
inline int vtkMIFloor(float x)
{
#if defined mips || defined sparc || defined __ppc__
  return (int)((unsigned int)(x + 2147483648.0) - 2147483648U);
#elif defined i386 || defined _M_IX86
  unsigned int hilo[2];
  *((double *)hilo) = x + 103079215104.0;  // (2**(52-16))*1.5
  return (int)((hilo[1]<<16)|(hilo[0]>>16));
#else
  return int(floor(x));
#endif
}


//----------------------------------------------------------------------------
// This templated function executes the filter for any type of data.
template <class T>
void vtkImageMutualInformationExecute(vtkImageMutualInformation *self,
				      T *in1Ptr, T *in2Ptr,
				      vtkImageData *inData1,
				      vtkImageData *inData2,
				      vtkImageData *outData,
				      int *outPtr,
				      double *NormalizedMI)
{
  int idX, idY, idZ;
  int iter, pmin0, pmax0, min0, max0, min1, max1, min2, max2;
  vtkIdType inInc0, inInc1, inInc2;
  T *tempPtr1, *tempPtr2;
  int *outPtrC;
  int  outIdx, outIdy, *outExtent;
  vtkIdType outIncs[3];
  vtkFloatingPointType *origin, *spacing;
  unsigned long count = 0;
  unsigned long target;
  int voxelCount = 0;

  int *xHist, *yHist; // arrays for the individual X and Y histograms

  vtkImageStencilData *stencil = self->GetStencil();

  // Zero count in every bin
  outData->GetExtent(min0, max0, min1, max1, min2, max2);
  memset((void *)outPtr, 0, 
         (max0-min0+1)*(max1-min1+1)*(max2-min2+1)*sizeof(int));
  // Set up and zero the image-specific arrays
  xHist = new int[max0-min0+1];
  yHist = new int[max1-min1+1];
  // Zero these arrays
  memset((void *)xHist, 0, (max0-min0+1)*sizeof(int));
  memset((void *)yHist, 0, (max1-min1+1)*sizeof(int));

    
  // Get information to march through data 
  inData1->GetWholeExtent(min0, max0, min1, max1, min2, max2);
  inData1->GetIncrements(inInc0, inInc1, inInc2);
  outExtent = outData->GetExtent();
  outData->GetIncrements(outIncs);
  origin = outData->GetOrigin();
  spacing = outData->GetSpacing();

  target = (unsigned long)((max2 - min2 + 1)*(max1 - min1 +1)/50.0);
  target++;

  // Loop through input pixels
  for (idZ = min2; idZ <= max2; idZ++)
    {
    for (idY = min1; idY <= max1; idY++)
      {
      if (!(count%target)) 
        {
        self->UpdateProgress(count/(50.0*target));
        }
      count++;

      // loop over stencil sub-extents
      iter = 0;
      if (self->GetReverseStencil())
        { // flag that we want the complementary extents
        iter = -1;
        }

      pmin0 = min0;
      pmax0 = max0;
      while ((stencil != 0 && 
              stencil->GetNextExtent(pmin0,pmax0,min0,max0,idY,idZ,iter)) ||
             (stencil == 0 && iter++ == 0))
        {
        // set up pointer for sub extent
        tempPtr1 = in1Ptr + (inInc2*(idZ - min2) +
			     inInc1*(idY - min1) +
			     (pmin0 - min0));
        tempPtr2 = in2Ptr + (inInc2*(idZ - min2) +
			     inInc1*(idY - min1) +
			     (pmin0 - min0));

        // accumulate over the sub extent
        for (idX = pmin0; idX <= pmax0; idX++)
          {
	  // compute the indices
	  outIdx = (int)*tempPtr1++;
	  outIdy = (int)*tempPtr2++;

	  // the imageA data are the 'x' dimension
	  // the imageB data are the 'y' dimension
	  // This is set up to leave 'top' and 'bottom' bins zeroed
	  if ((outIdx > outExtent[0]) && (outIdx < outExtent[1]) &&
	      (outIdy > outExtent[2]) && (outIdy < outExtent[3]))
	    {
            // In bin range
	    voxelCount++;
            outPtrC = outPtr + (outIdy * outIncs[1]) + outIdx;
	    ++(*outPtrC);
	    xHist[outIdx]++;
	    yHist[outIdy]++;
	    }
          }
        }
      }
    }
  
  double xEntropy, yEntropy, xyEntropy;
  xEntropy = yEntropy = xyEntropy = 0.0;

  for (outIdx = 0; outIdx < outExtent[1]; outIdx++)
    {
    if (xHist[outIdx]>0)
      xEntropy += (xHist[outIdx]*log((double)(xHist[outIdx])));
    }
  for (outIdy = 0; outIdy < outExtent[3]; outIdy++)
    {
    if (yHist[outIdy]>0)
      yEntropy += (yHist[outIdy]*log((double)(yHist[outIdy])));
    }

  outPtrC = outPtr;
  for (outIdx = 0; outIdx < outIncs[2]; outIdx++)
    {
    if ((*outPtrC)>0)
      xyEntropy += ((*outPtrC) * log((double)(*outPtrC)));
    outPtrC++;
    }

  xEntropy  = -   xEntropy / (double)voxelCount + log((double)voxelCount);
  yEntropy  = -   yEntropy / (double)voxelCount + log((double)voxelCount);
  xyEntropy = -  xyEntropy / (double)voxelCount + log((double)voxelCount);

  *NormalizedMI = (xEntropy+yEntropy)/xyEntropy;

  delete [] xHist;
  delete [] yHist;

}

        

//----------------------------------------------------------------------------
// This method is passed a input and output Data, and executes the filter
// algorithm to fill the output from the input.
// It just executes a switch statement to call the correct function for
// the Datas data types.
void vtkImageMutualInformation::ExecuteData(vtkDataObject *vtkNotUsed(out))
{
  void *inPtr1, *inPtr2;
  void *outPtr;

  vtkImageData *inData1 = this->GetInput1();
  if (inData1 == NULL)
    {
    vtkErrorMacro(<<"Input 0 must be specified.");
    return;
    }

  vtkImageData *inData2 = this->GetInput2();
  if (inData2 == NULL)
    {
    vtkErrorMacro(<<"Input 1 must be specified.");
    return;
    }

  vtkImageData *outData = this->GetOutput();
  
  vtkDebugMacro(<<"In vtkImageMutualInformation::ExecuteData");
  
  // We need to allocate our own scalars since we are overriding
  // the superclasses "Execute()" method.
  outData->SetExtent(outData->GetWholeExtent());
  outData->AllocateScalars();
  
  inPtr1 = inData1->GetScalarPointer();
  inPtr2 = inData2->GetScalarPointer();
  outPtr = outData->GetScalarPointer();
  
  // Components turned into x, y and z
  if (this->GetInput()->GetNumberOfScalarComponents() > 1)
    {
    vtkErrorMacro("This filter can only handle 1 scalar component.");
    return;
    }
  
  // this filter expects that output is type int.
  if (outData->GetScalarType() != VTK_INT)
    {
    vtkErrorMacro(<< "Execute: out ScalarType " << outData->GetScalarType()
                  << " must be int\n");
    return;
    }

  int wholeInExt1[6];
  int wholeInExt2[6];
  inData1->GetWholeExtent(wholeInExt1);
  inData2->GetWholeExtent(wholeInExt2);
  
  if ((wholeInExt1[0] - wholeInExt2[0]) | (wholeInExt1[1] - wholeInExt2[1]) |
      (wholeInExt1[2] - wholeInExt2[2]) | (wholeInExt1[3] - wholeInExt2[3]) |
      (wholeInExt1[4] - wholeInExt2[4]) | (wholeInExt1[5] - wholeInExt2[5]))
    {
    vtkErrorMacro(<<"Inputs must have matching wholeExtents.");
    return;
    }

  // this filter expects that inputs that have the same number of components
  if ((inData1->GetNumberOfScalarComponents() !=1) |
      (inData2->GetNumberOfScalarComponents() !=1))
    {
    vtkErrorMacro(<< "Execute: NumberOfScalarComponents must be 1 in all inputs");
    return;
    }

  // scalar types of the inputs must be the same
  if (inData1->GetScalarType() != inData2->GetScalarType())
    {
    vtkErrorMacro(<< "Execute: ScalarType Mismatch");
    return;
    }

  switch (inData1->GetScalarType())
    {
    vtkTemplateMacro8(vtkImageMutualInformationExecute, this, 
		      (VTK_TT *)(inPtr1), (VTK_TT *)(inPtr2),
		      inData1,inData2, 
		      outData, (int *)(outPtr),
		      &this->NormalizedMI);
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
      return;
    }
}


//----------------------------------------------------------------------------
void vtkImageMutualInformation::ExecuteInformation(vtkImageData **inputs, 
						   vtkImageData *output)
{
  // the two inputs are required to be of the same data type and extents.
  output->SetWholeExtent(this->ImageAComponentExtent[0],
			 this->ImageAComponentExtent[1],
			 this->ImageBComponentExtent[0],
			 this->ImageBComponentExtent[1],0,0);
  output->SetOrigin(this->ImageAComponentOrigin,
		    this->ImageBComponentOrigin,0.0f);
  output->SetSpacing(this->ImageAComponentSpacing,
		     this->ImageBComponentSpacing,1.0f);
  output->SetNumberOfScalarComponents(1);
  output->SetScalarType(VTK_INT);

  // need to set the spacing and origin of the stencil to match the output
  vtkImageStencilData *stencil = this->GetStencil();
  if (stencil)
    {
    stencil->SetSpacing(inputs[0]->GetSpacing());
    stencil->SetOrigin(inputs[0]->GetOrigin());
    }
}

//----------------------------------------------------------------------------
// Get ALL of the input.
void vtkImageMutualInformation::ComputeInputUpdateExtent(int inExt[6],
							 int outExt[6],
							 int vtkNotUsed(whichInput))
{
  int *wholeExtent = this->GetInput()->GetWholeExtent();

  memcpy(inExt, wholeExtent, 6*sizeof(int));
}


void vtkImageMutualInformation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Normalized MI: " << this->NormalizedMI << "\n";
  os << indent << "Stencil: " << this->GetStencil() << "\n";
  os << indent << "ReverseStencil: " << (this->ReverseStencil ?
                                         "On\n" : "Off\n");

  os << indent << "ComponentOrigin: ( "
     << this->ImageAComponentOrigin << ", "
     << this->ImageBComponentOrigin << ", "
     << 0 << " )\n";

  os << indent << "ComponentSpacing: ( "
     << this->ImageAComponentSpacing << ", "
     << this->ImageBComponentSpacing << ", "
     << 1 << " )\n";

  os << indent << "ComponentExtent: ( "
     << this->ImageAComponentExtent[0] << "," << this->ImageAComponentExtent[1] << " "
     << this->ImageBComponentExtent[0] << "," << this->ImageBComponentExtent[1] << " "
     << 0 << "," << 0 << " )\n";
}
