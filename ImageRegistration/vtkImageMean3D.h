/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImageMean3D.h,v $
  Language:  C++
  Date:      $Date: 2006/11/10 18:31:42 $
  Version:   $Revision: 1.2 $
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
// .NAME vtkImageMean3D - Mean Filter - this name should be changed!
// .SECTION Description

// vtkVectorSmooth a Smoothing filter used in the nonlinear warping system.
// 
// Replaces all vectors whose magnitude squared is greater than the threshold
// with a weighted smoothing of the surrounding vectors.
// Vectors shorter than the threshold are copied across.
// KernelSize may be specified.


#ifndef __vtkImageMean3D_h
#define __vtkImageMean3D_h

#include "vtkImageSpatialAlgorithm.h"

class VTK_EXPORT vtkImageMean3D : public vtkImageSpatialAlgorithm
{
public:
  static vtkImageMean3D *New();
  vtkTypeMacro(vtkImageMean3D,vtkImageSpatialAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // This method sets the size of the neighborhood.  It also sets the 
  // default middle of the neighborhood 
  void SetKernelSize(int size0, int size1, int size2);


  // Description:
  // This method sets the threshold above which smoothing is done.
  // For three component data, the square of the length of the vector is computed
  // and compared to this value to do the determination.
  vtkSetMacro(SmoothThreshold, float);
  vtkGetMacro(SmoothThreshold, float);

  // Description:
  // Sets the weighting to be applied to the center vector of the smoothing
  // neighborhood.
  vtkSetMacro(CenterWeighting, float);
  vtkGetMacro(CenterWeighting, float);

  // Description:
  // Sets the weighting to be applied to the surround of the smoothing
  // neighborhood.
  vtkSetMacro(SurroundWeighting, float);
  vtkGetMacro(SurroundWeighting, float);

  // Description:
  // Return the number of elements in the mean mask
  vtkGetMacro(NumberOfElements,int);
  
protected:
  vtkImageMean3D();
  ~vtkImageMean3D() {};
  vtkImageMean3D(const vtkImageMean3D&) {};
  void operator=(const vtkImageMean3D&) {};

  int NumberOfElements;

  float SmoothThreshold;
  float CenterWeighting;
  float SurroundWeighting;

  void ThreadedExecute(vtkImageData *inData, vtkImageData *outData, 
		       int extent[6], int id);

};

#endif



