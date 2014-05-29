/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImage3DNoiseSource.h,v $
  Language:  C++
  Date:      $Date: 2007/08/24 20:02:25 $
  Version:   $Revision: 1.3 $
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
// .NAME vtkImage3DNoiseSource - Create an image filled with noise.
// .SECTION Description
// vtkImage3DNoiseSource just produces images filled with noise.  The only
// option now is uniform noise specified by a min and a max.  There is one
// major problem with this source. Every time it executes, it will output
// different pixel values.  This has important implications when a stream
// requests overlapping regions.  The same pixels will have different values
// on different updates.


#ifndef __vtkImage3DNoiseSource_h
#define __vtkImage3DNoiseSource_h


#include "vtkImageSource.h"


class VTK_EXPORT vtkImage3DNoiseSource : public vtkImageSource 
{
public:
  vtkTypeMacro(vtkImage3DNoiseSource,vtkImageSource);
  static vtkImage3DNoiseSource *New();

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set/Get the minimum and maximum values for the generated noise.
  vtkSetMacro(Minimum, float);
  vtkGetMacro(Minimum, float);
  vtkSetMacro(Maximum, float);
  vtkGetMacro(Maximum, float);

  // Description:
  // Set how large of an image to generate.
  void SetWholeExtent(int xMinx, int xMax, int yMin, int yMax,
		      int zMin, int zMax);

  // Description:
  // Set the desired output scalar type.
  vtkSetMacro(OutputScalarType, int);
  vtkGetMacro(OutputScalarType, int);
  void SetOutputScalarTypeToDouble()
  {this->SetOutputScalarType(VTK_DOUBLE);}
  void SetOutputScalarTypeToFloat()
  {this->SetOutputScalarType(VTK_FLOAT);}
  void SetOutputScalarTypeToLong()
  {this->SetOutputScalarType(VTK_LONG);}
  void SetOutputScalarTypeToUnsignedLong()
  {this->SetOutputScalarType(VTK_UNSIGNED_LONG);};
  void SetOutputScalarTypeToInt()
  {this->SetOutputScalarType(VTK_INT);}
  void SetOutputScalarTypeToUnsignedInt()
  {this->SetOutputScalarType(VTK_UNSIGNED_INT);}
  void SetOutputScalarTypeToShort()
  {this->SetOutputScalarType(VTK_SHORT);}
  void SetOutputScalarTypeToUnsignedShort()
  {this->SetOutputScalarType(VTK_UNSIGNED_SHORT);}
  void SetOutputScalarTypeToChar()
  {this->SetOutputScalarType(VTK_CHAR);}
  void SetOutputScalarTypeToUnsignedChar()
  {this->SetOutputScalarType(VTK_UNSIGNED_CHAR);}

  // Description:
  // Set/Get the number of scalar components to generate.
  void SetNumberOfScalarComponents(int n);
  vtkGetMacro(NumberOfScalarComponents, int);

protected:
  vtkImage3DNoiseSource();
  ~vtkImage3DNoiseSource() {};

  virtual void ExecuteInformation();
  virtual void ExecuteData(vtkDataObject *data);

  float Minimum;
  float Maximum;
  int WholeExtent[6];
  int OutputScalarType;
  int NumberOfScalarComponents;

private:
  vtkImage3DNoiseSource(const vtkImage3DNoiseSource&); // Not implemented.
  void operator=(const vtkImage3DNoiseSource&); // Not implemented.
};


#endif

  
