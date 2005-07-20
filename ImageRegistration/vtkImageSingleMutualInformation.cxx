/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImageSingleMutualInformation.cxx,v $
  Language:  C++
  Date:      $Date: 2005/07/20 15:35:27 $
  Version:   $Revision: 1.3 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageSingleMutualInformation.h"

#include "vtkImageData.h"
#include "vtkImageStencilData.h"
#include "vtkObjectFactory.h"

#include <math.h>

vtkCxxRevisionMacro(vtkImageSingleMutualInformation, "$Revision: 1.3 $");
vtkStandardNewMacro(vtkImageSingleMutualInformation);

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageSingleMutualInformation::vtkImageSingleMutualInformation()
{
  this->ImageAComponentSpacing = 1.0;
  this->ImageAComponentOrigin = 0.0;
  this->ImageAComponentExtent[0] = 0;
  this->ImageAComponentExtent[1] = 255;
  
  this->ReverseStencil = 0;

  this->NormalizedMI = 0.0;
  this->MeanVoxel = 0.0;
}


//----------------------------------------------------------------------------
vtkImageSingleMutualInformation::~vtkImageSingleMutualInformation()
{
}

//----------------------------------------------------------------------------
void vtkImageSingleMutualInformation::SetImageAComponentExtent(int extent[2])
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
void vtkImageSingleMutualInformation::SetImageAComponentExtent(int min, int max)
{
  int extent[2];
  
  extent[0] = min;  extent[1] = max;
  this->SetImageAComponentExtent(extent);
}

//----------------------------------------------------------------------------
void vtkImageSingleMutualInformation::GetImageAComponentExtent(int extent[2])
{
  extent[0] = this->ImageAComponentExtent[0];
  extent[1] = this->ImageAComponentExtent[1];
}

//----------------------------------------------------------------------------
void vtkImageSingleMutualInformation::SetStencil(vtkImageStencilData *stencil)
{
  this->vtkProcessObject::SetNthInput(2, stencil); 
}


//----------------------------------------------------------------------------
vtkImageStencilData *vtkImageSingleMutualInformation::GetStencil()
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
void vtkImageSingleMutualInformationExecute(vtkImageSingleMutualInformation *self,
				      T *in1Ptr,
				      vtkImageData *inData1,
				      vtkImageData *outData,
				      int *outPtr,
				      double *NormalizedMI,
				      double *MeanVoxel)
{
  int idX, idY, idZ, idxC;
  int iter, pmin0, pmax0, min0, max0, min1, max1, min2, max2;
  int inInc0, inInc1, inInc2;
  T *tempPtr1;
  int *outPtrC;
  int  outIdx, outIdy, *outExtent, *outIncs;
  vtkFloatingPointType *origin, *spacing;
  float sum = 0.0;
  unsigned long totalVoxels = 0;
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
  outIncs = outData->GetIncrements();
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
        // accumulate over the sub extent
        for (idX = pmin0; idX <= pmax0; idX++)
          {
	  // compute the indices
	  sum += *tempPtr1;
	  totalVoxels++;
	  outIdx = (int)*tempPtr1++;
	  outIdy = 0;

	  // the imageA data are the 'x' dimension
	  // This is set up to leave 'top' and 'bottom' bins zeroed
	  if ((outIdx > outExtent[0]) && (outIdx < outExtent[1]))
	    {
            // In bin range
	    voxelCount++;
            outPtrC = outPtr + outIdx;
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

  xEntropy  = -   xEntropy * log((double)xEntropy)  / (double)voxelCount + log((double)voxelCount);
  yEntropy  = -   yEntropy * log((double)yEntropy)  / (double)voxelCount + log((double)voxelCount);
  xyEntropy = -  xyEntropy * log((double)xyEntropy) / (double)voxelCount + log((double)voxelCount);

  cout << "xEntropy "<<xEntropy<<" yEntropy " <<yEntropy<< " voxelCount "<<voxelCount<<"\n";
  cout << "totalVoxels "<<totalVoxels<< " xyEntropy "<<xyEntropy<<"\n";

  if(voxelCount > 0)
    {    
    *NormalizedMI = (xEntropy+yEntropy)/xyEntropy;
    *MeanVoxel = sum / voxelCount;
    }
  else
    {
    *NormalizedMI = 0;
    *MeanVoxel = 0;
    }

  delete [] xHist;
  delete [] yHist;

}

        

//----------------------------------------------------------------------------
// This method is passed a input and output Data, and executes the filter
// algorithm to fill the output from the input.
// It just executes a switch statement to call the correct function for
// the Datas data types.
void vtkImageSingleMutualInformation::ExecuteData(vtkDataObject *vtkNotUsed(out))
{
  void *inPtr1;
  void *outPtr;

  vtkImageData *inData1 = this->GetInput();
  if (inData1 == NULL)
    {
    vtkErrorMacro(<<"Input 0 must be specified.");
    return;
    }

  vtkImageData *outData = this->GetOutput();
  
  vtkDebugMacro(<<"In vtkImageSingleMutualInformation::ExecuteData");
  
  // We need to allocate our own scalars since we are overriding
  // the superclasses "Execute()" method.
  outData->SetExtent(outData->GetWholeExtent());
  outData->SetScalarType(VTK_INT);
  outData->AllocateScalars();
  
  inPtr1 = inData1->GetScalarPointer();
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
  inData1->GetWholeExtent(wholeInExt1);
  
  switch (inData1->GetScalarType())
    {
    vtkTemplateMacro7(vtkImageSingleMutualInformationExecute, this, 
		      (VTK_TT *)(inPtr1),
		      inData1, 
		      outData, (int *)(outPtr),
		      &this->NormalizedMI,
		      &this->MeanVoxel);
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
      return;
    }
}


//----------------------------------------------------------------------------
void vtkImageSingleMutualInformation::ExecuteInformation(vtkImageData **inputs, 
							 vtkImageData *output)
{

  // the two inputs are required to be of the same data type and extents.

  inputs[0]->Update();

  output->SetWholeExtent(this->ImageAComponentExtent[0],
			 this->ImageAComponentExtent[1],
			 this->ImageAComponentExtent[0],
			 this->ImageAComponentExtent[1],0,0);
  output->SetOrigin(this->ImageAComponentOrigin,
		    this->ImageAComponentOrigin,0.0f);
  output->SetSpacing(this->ImageAComponentSpacing,
		     this->ImageAComponentSpacing,0.0f);
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
void vtkImageSingleMutualInformation::ComputeInputUpdateExtent(int inExt[6],
							       int outExt[6],
							       int vtkNotUsed(whichInput))
{
  int *wholeExtent = this->GetInput()->GetWholeExtent();

  memcpy(inExt, wholeExtent, 6*sizeof(int));
}


void vtkImageSingleMutualInformation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Normalized MI: " << this->NormalizedMI << "\n";
  os << indent << "Mean Voxel: " << this->MeanVoxel << "\n";
  os << indent << "Stencil: " << this->GetStencil() << "\n";
  os << indent << "ReverseStencil: " << (this->ReverseStencil ?
                                         "On\n" : "Off\n");

  os << indent << "ComponentOrigin: ( "
     << this->ImageAComponentOrigin << ", "
     << 0 << " )\n";

  os << indent << "ComponentSpacing: ( "
     << this->ImageAComponentSpacing << ", "
     << 1 << " )\n";

  os << indent << "ComponentExtent: ( "
     << this->ImageAComponentExtent[0] << "," << this->ImageAComponentExtent[1] << " "
     << 0 << "," << 1 << " }\n";
}

