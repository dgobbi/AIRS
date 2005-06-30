/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImageCrossCorrelation.h,v $
  Language:  C++
  Date:      $Date: 2005/06/30 15:55:20 $
  Version:   $Revision: 1.2 $
  Thanks:    Thanks to Yves who developed this class.

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
// .NAME vtkImageCrossCorrelation - This computes the normalized cross-
//          correlation of two images.
// .SECTION Description
// vtkImageCrossCorrelation downsamples according to outputextent.
// KernelRadius is specified in "input extent" units.

#ifndef __vtkImageCrossCorrelation_h
#define __vtkImageCrossCorrelation_h

#ifndef vtkFloatingPointType

#define vtkFloatingPointType vtkFloatingPointType

typedef float vtkFloatingPointType;

#endif

#include "vtkImageTwoInputFilter.h"

class VTK_EXPORT vtkImageCrossCorrelation : public vtkImageTwoInputFilter
{
public:
  static vtkImageCrossCorrelation *New();
  vtkTypeMacro(vtkImageCrossCorrelation,vtkImageTwoInputFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set/Get the shrink factors
  vtkSetVector6Macro(OutputExtent, int);
  vtkGetVector6Macro(OutputExtent, int);


  // Description:
  // Set/Get the kernel size
  vtkSetVector3Macro(KernelRadius, float);
  vtkGetVector3Macro(KernelRadius, float);

  vtkGetVector3Macro(ShrinkFactors, float);
    
protected:
  vtkImageCrossCorrelation();
  ~vtkImageCrossCorrelation() {};
  vtkImageCrossCorrelation(const vtkImageCrossCorrelation&) {};
  void operator=(const vtkImageCrossCorrelation&) {};

  float ShrinkFactors[3];
  float KernelRadius[3];
  int OutputExtent[6];

  void ExecuteInformation(vtkImageData **inDatas, vtkImageData *outData);
  void ComputeInputUpdateExtent(int inExt[6], int outExt[6], 
				int vtkNotUsed(whichInput));
  void ExecuteInformation(){this->vtkImageTwoInputFilter::ExecuteInformation();};
  void ThreadedExecute(vtkImageData **inDatas, vtkImageData *outData,
		       int extent[6], int id);

};

#endif



