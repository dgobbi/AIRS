/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImageFloodFill.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkImageFloodFill.h"

//#include "vtkImageData.h"
#include "vtkImageProgressIterator.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkImageStencilData.h"

#include "vtkImageImport.h"
#include "vtkImageShiftScale.h"

#include <stack>
#include "math.h"

vtkCxxRevisionMacro(vtkImageFloodFill, "$Revision: 1.4 $");
vtkStandardNewMacro(vtkImageFloodFill);
vtkCxxSetObjectMacro(vtkImageFloodFill, SeedPoints, vtkPoints);

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageFloodFill::vtkImageFloodFill()
{
  this->UpperThreshold = VTK_LARGE_FLOAT;
  this->LowerThreshold = -VTK_LARGE_FLOAT;
  this->SeedPoints = 0;
  this->ReplaceIn = 0;
  this->InValue = 0.0;
  this->ReplaceOut = 0;
  this->OutValue = 0.0;

  this->FloodExtent[0] = -VTK_INT_MAX;
  this->FloodExtent[2] = -VTK_INT_MAX;
  this->FloodExtent[4] = -VTK_INT_MAX;

  this->FloodExtent[1] = VTK_INT_MAX;
  this->FloodExtent[3] = VTK_INT_MAX;
  this->FloodExtent[5] = VTK_INT_MAX;

  this->FloodBounds[0] = VTK_INT_MAX;
  this->FloodBounds[2] = VTK_INT_MAX;
  this->FloodBounds[4] = VTK_INT_MAX;

  this->FloodBounds[1] = -VTK_INT_MAX;
  this->FloodBounds[3] = -VTK_INT_MAX;
  this->FloodBounds[5] = -VTK_INT_MAX;

  this->ReverseStencil = 0;
  this->ActiveComponent = 0;

  this->ImageMask = vtkImageData::New();
  this->ImageMask->SetScalarTypeToUnsignedChar();

  this->NumberOfInVoxels = 0;
}

//----------------------------------------------------------------------------
vtkImageFloodFill::~vtkImageFloodFill()
{
  if (this->SeedPoints)
    {
    this->SeedPoints->Delete();
    }
  this->ImageMask->Delete();
}


//----------------------------------------------------------------------------
void vtkImageFloodFill::SetInValue(double val)
{
  if (val != this->InValue || this->ReplaceIn != 1)
    {
    this->InValue = val;
    this->ReplaceIn = 1;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkImageFloodFill::SetOutValue(double val)
{
  if (val != this->OutValue || this->ReplaceOut != 1)
    {
    this->OutValue = val;
    this->ReplaceOut = 1;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
// The values greater than or equal to the value match.
void vtkImageFloodFill::ThresholdByUpper(double thresh)
{
  if (this->LowerThreshold != thresh || this->UpperThreshold < VTK_LARGE_FLOAT)
    {
    this->LowerThreshold = thresh;
    this->UpperThreshold = VTK_LARGE_FLOAT;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
// The values less than or equal to the value match.
void vtkImageFloodFill::ThresholdByLower(double thresh)
{
  if (this->UpperThreshold != thresh || this->LowerThreshold > -VTK_LARGE_FLOAT)
    {
    this->UpperThreshold = thresh;
    this->LowerThreshold = -VTK_LARGE_FLOAT;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
// The values in a range (inclusive) match
void vtkImageFloodFill::ThresholdBetween(double lower, double upper)
{
  if (this->LowerThreshold != lower || this->UpperThreshold != upper)
    {
    this->LowerThreshold = lower;
    this->UpperThreshold = upper;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkImageFloodFill::SetStencil(vtkImageStencilData *stencil)
{
  this->vtkProcessObject::SetNthInput(1, stencil); 
}

//----------------------------------------------------------------------------
vtkImageStencilData *vtkImageFloodFill::GetStencil()
{
  if (this->NumberOfInputs < 2) 
    { 
    return NULL;
    }
  else
    {
    return (vtkImageStencilData *)(this->Inputs[1]); 
    }
}

//----------------------------------------------------------------------------
unsigned long vtkImageFloodFill::GetMTime() 
{
  unsigned long mTime = this->MTime.GetMTime();
  unsigned long pointsMTime;

  if (this->SeedPoints)
    {
    pointsMTime = this->SeedPoints->GetMTime();
    mTime = ( pointsMTime > mTime ? pointsMTime : mTime );
    }

  return mTime;
}

//----------------------------------------------------------------------------
void vtkImageFloodFill::ComputeInputUpdateExtent(int inExt[6], 
						 int outExt[6])
{
  int extent[6];
  this->GetFloodExtent(extent);
  this->GetInput()->GetWholeExtent(inExt);
  int i;

  // Expand the extent to match the outExt, but only if 
  // ReplaceOut is set
  if (this->ReplaceOut)
    {
    for (i = 0; i < 3; i++)
      {
      if (extent[2*i] > outExt[2*i])
	{
	extent[2*i] = outExt[2*i];
	}
      if (extent[2*i+1] < outExt[2*i+1])
	{
	extent[2*i+1] = outExt[2*i+1];
	}
      }
    }

  // Clip the inExt to the bounds of the extent
  for (i = 0; i < 3; i++)
    {
    if (extent[2*i] > inExt[2*i+1] || extent[2*i+1] < inExt[2*i])
      { // extents don't intersect
      inExt[2*i] = inExt[2*i+1]+1; // signal for empty extent
      continue;
      }
    if (inExt[2*i] < extent[2*i])
      {
      inExt[2*i] = extent[2*i];
      }
    if (inExt[2*i+1] > extent[2*i+1])
      {
      inExt[2*i+1] = extent[2*i+1];
      }
    }
}

//----------------------------------------------------------------------------
void vtkImageFloodFill::ExecuteInformation(vtkImageData *inData, 
                                           vtkImageData *outData)
{
  // need to set the spacing and origin of the stencil to match the output
  vtkImageStencilData *stencil = this->GetStencil();
  if (stencil)
    {
    stencil->SetSpacing(inData->GetSpacing());
    stencil->SetOrigin(inData->GetOrigin());
    }
}

//----------------------------------------------------------------------------
// seed struct: just a set of indices
class vtkFloodFillSeed
{
public:
  vtkFloodFillSeed() {
    store[0]=0; store[1]=0; store[2]=0; }; 
  vtkFloodFillSeed(int i, int j, int k) { 
    store[0]=i; store[1]=j; store[2]=k; };
  vtkFloodFillSeed(const vtkFloodFillSeed &seed) { 
    store[0]=seed.store[0]; store[1]=seed.store[1]; store[2]=seed.store[2]; };
  const int &operator[](int i) const { return store[i]; };
  const vtkFloodFillSeed &operator=(const vtkFloodFillSeed seed) {
    store[0]=seed.store[0]; store[1]=seed.store[1]; store[2]=seed.store[2]; 
    return *this; };

private:    
  int store[3];
};

//----------------------------------------------------------------------------
template <class IT>
void vtkImageFloodFillThresholds(vtkImageFloodFill *self,
                                 vtkImageData *inData,
                                 IT &lowerThreshold,
                                 IT &upperThreshold)
{
  // Make sure the thresholds are valid for the input scalar range
  if (self->GetLowerThreshold() < inData->GetScalarTypeMin())
    {
    lowerThreshold = (IT)inData->GetScalarTypeMin();
    }
  else 
    {
    if (self->GetLowerThreshold() > inData->GetScalarTypeMax())
      {
      lowerThreshold = (IT)inData->GetScalarTypeMax();
      }
    else
      {
      lowerThreshold = (IT)self->GetLowerThreshold();
      }
    }
  if (self->GetUpperThreshold() > inData->GetScalarTypeMax())
    {
    upperThreshold = (IT)inData->GetScalarTypeMax();
    }
  else 
    {
    if (self->GetUpperThreshold() < inData->GetScalarTypeMin())
      {
      upperThreshold = (IT)inData->GetScalarTypeMin();
      }
    else
      {
      upperThreshold = (IT)self->GetUpperThreshold();
      }
    }
}

//----------------------------------------------------------------------------
template <class OT>
void vtkImageFloodFillValues(vtkImageFloodFill *self,
                             vtkImageData *outData,
                             OT &inValue,
                             OT &outValue)
{
  // Make sure the replacement values are within the output scalar range
  if (self->GetInValue() < outData->GetScalarTypeMin())
    {
    inValue = (OT) outData->GetScalarTypeMin();
    }
  else 
    {
    if (self->GetInValue() > outData->GetScalarTypeMax())
      {
      inValue = (OT) outData->GetScalarTypeMax();
      }
    else
      {
      inValue = (OT) self->GetInValue();
      }
    }
  if (self->GetOutValue() > outData->GetScalarTypeMax())
    {
    outValue = (OT) outData->GetScalarTypeMax();
    }
  else 
    {
    if (self->GetOutValue() < outData->GetScalarTypeMin())
      {
      outValue = (OT) outData->GetScalarTypeMin();
      }
    else
      {
      outValue = (OT) self->GetOutValue();
      }
    }
}

//----------------------------------------------------------------------------
void vtkImageFloodFillApplyStencil(vtkImageFloodFill *self,
                                   unsigned char *ptr,
                                   vtkImageStencilData *stencil,
                                   int extent[6])
{
  // if there is a stencil, use it to initialize the mask
  int idX, idY, idZ;

  for (idZ = extent[4]; idZ <= extent[5]; idZ++)
    {
    for (idY = extent[2]; idY <= extent[3]; idY++)
      {
      int iter = (self->GetReverseStencil() ? -1 : 0);
      int c1 = extent[0];
      int r1 = extent[0];
      int r2 = extent[1];
      for (;;)
        {
        if (stencil->GetNextExtent(r1, r2, extent[0], extent[1],
                                   idY, idZ, iter))
          {
          for (idX = c1; idX < r1; idX++)
            {
            *ptr++ = 255;
            }
          for (idX = r1; idX <= r2; idX++)
            {
            *ptr++ = 0;
            }
          c1 = r2 + 1;
          }
        else
          {
          for (idX = c1; idX <= extent[1]; idX++)
            {
            *ptr++ = 1;
            }
          break;
          }
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkImageFloodFillApplyMaskStencil(vtkImageFloodFill *self,
				       unsigned char *ptr,
				       vtkImageStencilData *stencil,
				       int extent[6])
{
  // if there is a stencil, use it to initialize the mask
  int idX, idY, idZ;
  unsigned char value;

  for (idZ = extent[4]; idZ <= extent[5]; idZ++)
    {
    for (idY = extent[2]; idY <= extent[3]; idY++)
      {
      int iter = (self->GetReverseStencil() ? -1 : 0);
      int c1 = extent[0];
      int r1 = extent[0];
      int r2 = extent[1];
      for (;;)
        {
        if (stencil->GetNextExtent(r1, r2, extent[0], extent[1],
                                   idY, idZ, iter))
          {
          for (idX = c1; idX < r1; idX++)
            {
            *ptr++ = 0;
            }
          for (idX = r1; idX <= r2; idX++)
            {
            //*ptr++ = *ptr;
	    value = *ptr;
	    *ptr++ = value;
            }
          c1 = r2 + 1;
          }
        else
          {
          for (idX = c1; idX <= extent[1]; idX++)
            {
            *ptr++ = 0;
            }
          break;
          }
        }
      }
    }
}


//----------------------------------------------------------------------------
// This templated function executes the filter for any type of data.
template <class IT, class OT>
void vtkImageFloodFillExecute(vtkImageFloodFill *self,
                              vtkImageData *inData,
                              vtkImageData *outData, 
			      vtkImageData *maskData,
                              int outExt[6], int id,
			      IT *inPtr, OT *outPtr)
{
  vtkImageIterator<IT> inIt(inData, outExt);
  vtkImageProgressIterator<OT> outIt(outData, outExt, self, id);
  IT  lowerThreshold;
  IT  upperThreshold;
  int replaceIn = self->GetReplaceIn();
  OT  inValue;
  int replaceOut = self->GetReplaceOut();
  OT  outValue;
  int nComponents = outData->GetNumberOfScalarComponents();
  int activeComponent = self->GetActiveComponent();
  activeComponent = activeComponent % nComponents;
  while (activeComponent < 0)
    {
    activeComponent += nComponents;
    }
  
  vtkImageFloodFillThresholds(self,inData,lowerThreshold,upperThreshold);
  vtkImageFloodFillValues(self,outData,inValue,outValue);
  
  // Set the "outside" with either the input or the outValue
  while (!outIt.IsAtEnd())
    {
    IT* inSI = inIt.BeginSpan();
    OT* outSI = outIt.BeginSpan();
    OT* outSIEnd = outIt.EndSpan();

    if (replaceOut)
      {
      if (nComponents == 1)
	{
	while (outSI < outSIEnd)
	  {
	  *outSI++ = outValue;
	  }
	}
      else
	{
	// only color the active component, copy the rest
	while (outSI < outSIEnd)
	  {
	  int jj = 0;
	  while (jj < activeComponent)
	    {
	    *outSI++ = (OT)(*inSI++);
	    jj++;
	    }
	  *outSI++ = outValue;
	  inSI++;
	  jj++;
	  while (jj < nComponents)
	    {
	    *outSI++ = (OT)(*inSI++);
	    jj++;
	    }
	  }
	}
      }
    else
      {
      while (outSI < outSIEnd)
	{
        *outSI++ = (OT)(*inSI++);
	}
      }
    inIt.NextSpan();
    outIt.NextSpan();
    }

  // Get the extent for the flood fill, and clip with the input extent
  int extent[6];
  int inExt[6];
  self->GetFloodExtent(extent);
  inData->GetExtent(inExt);
  int outCheck = 0;
  for (int ii = 0; ii < 3; ii++)
    {
    if (extent[2*ii] > inExt[2*ii+1] || extent[2*ii+1] < inExt[2*ii])
      { // extents don't intersect, we're done
      return;
      }
    if (extent[2*ii] < inExt[2*ii])
      {
      extent[2*ii] = inExt[2*ii];
      }
    if (extent[2*ii+1] > inExt[2*ii+1])
      {
      extent[2*ii+1] = inExt[2*ii+1];
      }
    // check against output extent
    if (extent[2*ii] < outExt[2*ii] || extent[2*ii+1] > outExt[2*ii+1])
      {
      outCheck = 1;
      }
    }

  // I think these are not the pointers we want
  inPtr = (IT *)inData->GetScalarPointerForExtent(extent);
  outPtr = (OT *)outData->GetScalarPointerForExtent(extent);
  int fullsize = ((vtkIdType)(extent[1]-extent[0]+1)*
		  (vtkIdType)(extent[3]-extent[2]+1)*
		  (vtkIdType)(extent[5]-extent[4]+1));


  // Setup the mask
  maskData->SetWholeExtent(inData->GetWholeExtent());
  maskData->SetExtent(inData->GetWholeExtent());
  maskData->SetOrigin(inData->GetOrigin());
  maskData->SetSpacing(inData->GetSpacing());

  unsigned char *maskPtr;
  maskPtr = (unsigned char *)maskData->GetScalarPointerForExtent(outExt);

  // Adjust pointers to active component
  inPtr += activeComponent;
  outPtr += activeComponent;

  if (self->GetStencil() == 0)
    {
    memset(maskPtr,0,fullsize);
    }
  else
    {
    vtkImageFloodFillApplyStencil(self,maskPtr,self->GetStencil(),extent);
    }
  
  // Perform the flood fill within the extent
  vtkIdType inInc[3];
  vtkIdType outInc[3];
  vtkIdType maskInc[3];
  inData->GetIncrements(inInc);
  outData->GetIncrements(outInc);
  maskInc[0] = 1;
  maskInc[1] = (extent[1]-extent[0]+1);
  maskInc[2] = maskInc[1]*(extent[3]-extent[2]+1);

#if (VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION <= 3) 
  float spacing[3];
  float origin[3];
#else
  double spacing[3];
  double origin[3];
#endif
  outData->GetSpacing(spacing);
  outData->GetOrigin(origin);

  // create the seed stack
  // stack has methods empty(), top(), push(), and pop()
  std::stack<vtkFloodFillSeed> seedStack;

  // initialize with the seeds provided by the user
  vtkPoints *points = self->GetSeedPoints();
  if (points == 0)
    { // no seeds!
    delete [] maskPtr;
    return;
    }
  double point[3];
  int nPoints = points->GetNumberOfPoints();
  for (int p = 0; p < nPoints; p++)
    {
    points->GetPoint(p,point);
    vtkFloodFillSeed seed = vtkFloodFillSeed(
		   int(floor((point[0]-origin[0])/spacing[0]+0.5))-extent[0],
		   int(floor((point[1]-origin[1])/spacing[1]+0.5))-extent[2],
		   int(floor((point[2]-origin[2])/spacing[2]+0.5))-extent[4]);

    if (seed[0] >= 0 && seed[0] <= extent[1]-extent[0] &&
	seed[1] >= 0 && seed[1] <= extent[3]-extent[2] &&
	seed[2] >= 0 && seed[2] <= extent[5]-extent[4])
      {
      seedStack.push(seed);
      }
    }

  int counter = 0;
  int FloodBounds[6];
  if (!seedStack.empty())
    {
    counter += 1; // count the seed
    vtkFloodFillSeed seed = seedStack.top();
    FloodBounds[0] = seed[0];
    FloodBounds[1] = seed[0];
    FloodBounds[2] = seed[1];
    FloodBounds[3] = seed[1];
    FloodBounds[4] = seed[2];
    FloodBounds[5] = seed[2];
    }
  else
    {
    FloodBounds[0] = VTK_INT_MAX;
    FloodBounds[2] = VTK_INT_MAX;
    FloodBounds[4] = VTK_INT_MAX;
    FloodBounds[1] = -VTK_INT_MAX;
    FloodBounds[3] = -VTK_INT_MAX;
    FloodBounds[5] = -VTK_INT_MAX;
    }
  
  while (!seedStack.empty())
    {
    vtkFloodFillSeed seed = seedStack.top();
    seedStack.pop();

    unsigned char *maskPtr1 = maskPtr + (seed[0]*maskInc[0] + 
					 seed[1]*maskInc[1] + 
					 seed[2]*maskInc[2]);

    if (*maskPtr1)
      {
      continue;
      }
    *maskPtr1 = 255;

    IT *inPtr1 = inPtr + ((vtkIdType)seed[0]*inInc[0] + 
			  (vtkIdType)seed[1]*inInc[1] + 
			  (vtkIdType)seed[2]*inInc[2]);
    IT temp = *inPtr1;

    if (lowerThreshold <= temp && temp <= upperThreshold)
      {
      // match
      OT *outPtr1 = outPtr + ((vtkIdType)seed[0]*outInc[0] + 
                              (vtkIdType)seed[1]*outInc[1] + 
                              (vtkIdType)seed[2]*outInc[2]);

      if (outCheck == 0 ||
	  seed[0] >= 0 && seed[0] <= outExt[1]-outExt[0] &&
          seed[1] >= 0 && seed[1] <= outExt[3]-outExt[2] &&
          seed[2] >= 0 && seed[2] <= outExt[5]-outExt[4])
        {
        if (replaceIn)
          {
          *outPtr1 = inValue;
          }
        else
          {
          *outPtr1 = (OT)(temp);
          }
        }
      // add new seeds to the stack
      counter += 1;

      // keep track of the flood fill bounds
      // X
      //FloodBounds[0] = FloodBounds[0] < seed[0] ? (FloodBounds[0] : seed[0]);

      if (seed[0] < FloodBounds[0])
	{
	FloodBounds[0] = seed[0];
	}
      else if (seed[0] > FloodBounds[1])
	{
	FloodBounds[1] = seed[0];
	}
      // Y
      if (seed[1] < FloodBounds[2])
	{
	FloodBounds[2] = seed[1];
	}
      else if (seed[1] > FloodBounds[3])
	{
	FloodBounds[3] = seed[1];
	}
      // Z
      if (seed[2] < FloodBounds[4])
	{
	FloodBounds[4] = seed[2];
	}
      else if (seed[2] > FloodBounds[5])
	{
	FloodBounds[5] = seed[2];
	}

      // push the seend
      if (seed[2] > 0 && *(maskPtr1 - maskInc[2]) == 0)
	{
	seedStack.push(vtkFloodFillSeed(seed[0],seed[1],seed[2]-1));
	}
      if (seed[2] < extent[5]-extent[4] && *(maskPtr1 + maskInc[2]) == 0)
	{
	seedStack.push(vtkFloodFillSeed(seed[0],seed[1],seed[2]+1));
	}
      if (seed[1] > 0 && *(maskPtr1 - maskInc[1]) == 0)
	{
	seedStack.push(vtkFloodFillSeed(seed[0],seed[1]-1,seed[2]));
	}
      if (seed[1] < extent[3]-extent[2] && *(maskPtr1 + maskInc[1]) == 0)
	{
	seedStack.push(vtkFloodFillSeed(seed[0],seed[1]+1,seed[2]));
	}
      if (seed[0] > 0 && *(maskPtr1 - maskInc[0]) == 0)
	{
	seedStack.push(vtkFloodFillSeed(seed[0]-1,seed[1],seed[2]));
	}
      if (seed[0] < extent[1]-extent[0] && *(maskPtr1 + maskInc[0]) == 0)
	{
	seedStack.push(vtkFloodFillSeed(seed[0]+1,seed[1],seed[2]));
	}
      }
    }

  self->SetNumberOfInVoxels(counter);
  self->SetFloodBounds(FloodBounds);

  if (self->GetStencil())
    {
    vtkImageFloodFillApplyMaskStencil(self,maskPtr,self->GetStencil(),extent);
    }
}

//----------------------------------------------------------------------------
// This is the superclasses style of Execute method.  Convert it into
// an imaging style Execute method, and don't multi-thread.
void vtkImageFloodFill::ExecuteData(vtkDataObject *out)
{
  vtkImageData *outData = this->AllocateOutputData(out);
  vtkImageData *inData = this->GetInput();
  vtkImageData *maskData = this->GetImageMask();
  int outExt[6];
  outData->GetUpdateExtent(outExt);
  void *inPtr = inData->GetScalarPointerForExtent(outExt);
  void *outPtr = outData->GetScalarPointerForExtent(outExt);
  int id = 0; // not multi-threaded

  if (inData->GetScalarType() != outData->GetScalarType())
    { 
    vtkErrorMacro("Execute: Output ScalarType " 
                  << outData->GetScalarType()
                  << ", must Input ScalarType "
                  << inData->GetScalarType());
    return;
    }

  switch (inData->GetScalarType())
    {
#if (VTK_MAJOR_VERSION < 5)
    vtkTemplateMacro8(vtkImageFloodFillExecute, this, inData,
                      outData, maskData, outExt, id, 
                      static_cast<VTK_TT *>(inPtr),
                      static_cast<VTK_TT *>(outPtr));
#else
    vtkTemplateMacro(
      vtkImageFloodFillExecute(this, inData,
			       outData, maskData, outExt, id, 
			       static_cast<VTK_TT *>(inPtr),
			       static_cast<VTK_TT *>(outPtr)));
#endif
    default:
      vtkErrorMacro(<< "Execute: Unknown input ScalarType");
      return;
    }
}

//----------------------------------------------------------------------------
void vtkImageFloodFill::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "InValue: " << this->InValue << "\n";
  os << indent << "OutValue: " << this->OutValue << "\n";
  os << indent << "LowerThreshold: " << this->LowerThreshold << "\n";
  os << indent << "UpperThreshold: " << this->UpperThreshold << "\n";
  os << indent << "ReplaceIn: " << this->ReplaceIn << "\n";
  os << indent << "ReplaceOut: " << this->ReplaceOut << "\n";
  os << indent << "FloodExtent: " << this->FloodExtent[0] << " " <<
    this->FloodExtent[1] << " " << this->FloodExtent[2] << " " <<
    this->FloodExtent[3] << " " << this->FloodExtent[4] << " " <<
    this->FloodExtent[5] << "\n";
  os << indent << "SeedPoints: " << this->SeedPoints << "\n";
  if (this->SeedPoints) 
    {
    this->SeedPoints->PrintSelf(os,indent.GetNextIndent());
    }
  os << indent << "Stencil: " << this->GetStencil() << "\n";
  os << indent << "ReverseStencil: " << (this->ReverseStencil ?
                                         "On\n" : "Off\n");
  os << indent << "ActiveComponent: " << this->ActiveComponent << "\n";
}

