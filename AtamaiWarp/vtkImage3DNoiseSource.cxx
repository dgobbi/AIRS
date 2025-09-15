/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImage3DNoiseSource.cxx,v $
  Language:  C++
  Date:      $Date: 2007/08/24 20:02:25 $
  Version:   $Revision: 1.4 $
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
#include <stdlib.h>
#include "vtkMath.h"
#include "vtkImageData.h"
#include "vtkImage3DNoiseSource.h"
#include "vtkObjectFactory.h"

//----------------------------------------------------------------------------
vtkImage3DNoiseSource* vtkImage3DNoiseSource::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkImage3DNoiseSource");
  if(ret)
  {
    return (vtkImage3DNoiseSource*)ret;
  }
  // If the factory was unable to create the object, then create it here.
  return new vtkImage3DNoiseSource;
}

//----------------------------------------------------------------------------
vtkImage3DNoiseSource::vtkImage3DNoiseSource()
{
  this->Minimum = 0.0;
  this->Maximum = 10.0;
  this->OutputScalarType = -1;
  this->WholeExtent[0] = 0;  this->WholeExtent[1] = 255;
  this->WholeExtent[2] = 0;  this->WholeExtent[3] = 255;
  this->WholeExtent[4] = 0;  this->WholeExtent[5] = 0;
  this->NumberOfScalarComponents = 1;
}

//----------------------------------------------------------------------------
void vtkImage3DNoiseSource::SetWholeExtent(int xMin, int xMax,
                                           int yMin, int yMax,
                                           int zMin, int zMax)
{
  int modified = 0;

  if (this->WholeExtent[0] != xMin)
  {
    modified = 1;
    this->WholeExtent[0] = xMin ;
  }
  if (this->WholeExtent[1] != xMax)
  {
    modified = 1;
    this->WholeExtent[1] = xMax ;
  }
  if (this->WholeExtent[2] != yMin)
  {
    modified = 1;
    this->WholeExtent[2] = yMin ;
  }
  if (this->WholeExtent[3] != yMax)
  {
    modified = 1;
    this->WholeExtent[3] = yMax ;
  }
  if (this->WholeExtent[4] != zMin)
  {
    modified = 1;
    this->WholeExtent[4] = zMin ;
  }
  if (this->WholeExtent[5] != zMax)
  {
    modified = 1;
    this->WholeExtent[5] = zMax ;
  }
  if (modified)
  {
    this->Modified();
  }
}

//----------------------------------------------------------------------------
void vtkImage3DNoiseSource::SetNumberOfScalarComponents(int num)
{
  this->NumberOfScalarComponents = num;
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkImage3DNoiseSource::ExecuteInformation()
{
  vtkImageData *output = this->GetOutput();

  output->SetWholeExtent(this->WholeExtent);
  output->SetScalarType(this->OutputScalarType);
  output->SetNumberOfScalarComponents(this->NumberOfScalarComponents);
}

//----------------------------------------------------------------------------
template <class T>
void vtkImage3DNoiseExecute(vtkImage3DNoiseSource *self, vtkImageData *data, T *outPtr)
{
  int idxR, idxY, idxZ;
  int maxY, maxZ;
  vtkIdType outIncX, outIncY, outIncZ;
  int rowLength;
  int *outExt;
  unsigned long count = 0;
  unsigned long target;
  double min = self->GetMinimum();
  double max = self->GetMaximum();

  outExt = data->GetExtent();

  // find the region to loop over
  rowLength = (outExt[1] - outExt[0]+1)*self->GetNumberOfScalarComponents();
  maxY = outExt[3] - outExt[2];
  maxZ = outExt[5] - outExt[4];

  // Get increments to march through data
  data->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);
  outPtr = (T *) data->GetScalarPointer(outExt[0],outExt[2],outExt[4]);

  target = (unsigned long)((maxZ+1)*(maxY+1)/50.0);
  target++;

  // Loop through output pixels
  for (idxZ = 0; idxZ <= maxZ; idxZ++)
  {
    for (idxY = 0; !self->AbortExecute && idxY <= maxY; idxY++)
    {
      if (!(count%target))
      {
	self->UpdateProgress(count/(50.0*target));
      }
      count++;
      for (idxR = 0; idxR < rowLength; idxR++)
      {
	// Pixel operation
	*outPtr = (T) (min + (max - min) * vtkMath::Random());
	outPtr++;
      }
      outPtr += outIncY;
    }
    outPtr += outIncZ;
  }
}

//----------------------------------------------------------------------------
void vtkImage3DNoiseSource::ExecuteData(vtkDataObject *output)
{
  vtkImageData *data = this->AllocateOutputData(output);
  void *outPtr = data->GetScalarPointer();

  switch (this->OutputScalarType)
  {
#if (VTK_MAJOR_VERSION < 5)
      vtkTemplateMacro3(vtkImage3DNoiseExecute, this, data, (VTK_TT *)(outPtr));
#else
      vtkTemplateMacro(vtkImage3DNoiseExecute(this, data, (VTK_TT *)(outPtr)));
#endif
    default:
      vtkErrorMacro(<< "Execute: Unknown output ScalarType");
  }
}

//----------------------------------------------------------------------------
void vtkImage3DNoiseSource::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkImageSource::PrintSelf(os,indent);

  os << indent << "Minimum: " << this->Minimum << "\n";
  os << indent << "Maximum: " << this->Maximum << "\n";
  os << indent << "Output Scalar Type: " << this->OutputScalarType << "\n";
  os << indent << "NumberOfScalarComponents: " <<
    this->NumberOfScalarComponents << "\n";
}

