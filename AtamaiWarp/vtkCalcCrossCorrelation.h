/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCalcCrossCorrelation.h,v $
  Language:  C++
  Date:      $Date: 2007/08/24 20:02:25 $
  Version:   $Revision: 1.3 $


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
// .NAME vtkCalcCrossCorrelation - normalized cross correlation of volumes
// .SECTION Description
// vtkCalcCrossCorrelation - normalized cross correlation

// .SECTION Caveats
// Volumes must be same type, extent, spacing


#ifndef __vtkCalcCrossCorrelation_h
#define __vtkCalcCrossCorrelation_h

#include "vtkProcessObject.h"

class vtkImageData;
class vtkImageStencilData;
class vtkInformation;

class VTK_EXPORT vtkCalcCrossCorrelation : public vtkProcessObject
{
public:
  vtkTypeMacro(vtkCalcCrossCorrelation,vtkProcessObject);
  static vtkCalcCrossCorrelation *New();

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Use a stencil to specify which voxels to accumulate.
  void SetStencil(vtkImageStencilData *stencil);
  vtkImageStencilData *GetStencil();

  // Description:
  // Reverse the stencil.
  vtkSetMacro(ReverseStencil, int);
  vtkBooleanMacro(ReverseStencil, int);
  vtkGetMacro(ReverseStencil, int);

  // Description:
  // Compute and return the cross correlation.
  double GetCrossCorrelation() {this->Update(); return this->CrossCorrelation;}

  void Update();
  
  void SetInput1(vtkImageData *input1);
  void SetInput2(vtkImageData *input2);
  vtkImageData *GetInput1();
  vtkImageData *GetInput2();

protected:
  vtkCalcCrossCorrelation();
  ~vtkCalcCrossCorrelation();
 
  virtual int FillInputPortInformation(int, vtkInformation*);

  void Execute();

  int ReverseStencil;

  double  CrossCorrelation;
  vtkTimeStamp ExecuteTime;

private:
  vtkCalcCrossCorrelation(const vtkCalcCrossCorrelation&); // Not implemented.
  void operator=(const vtkCalcCrossCorrelation&); // Not implemented.

};

#endif








