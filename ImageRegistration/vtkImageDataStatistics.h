/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImageDataStatistics.h,v $
  Language:  C++
  Date:      $Date: 2004/07/13 14:43:11 $
  Version:   $Revision: 1.1 $
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
// .NAME vtkImageDataStatistics - Computes displacement grid statistics.
//
// .SECTION Description
// Takes a 3 scalar component vtkImageData as input.
// - computes mean and standard deviation of vector magnitude.

#ifndef __vtkImageDataStatistics_h
#define __vtkImageDataStatistics_h

#include "vtkProcessObject.h"
#include "vtkImageData.h"

class VTK_EXPORT vtkImageDataStatistics : public vtkProcessObject
{
public:
  static vtkImageDataStatistics *New();

  vtkTypeMacro(vtkImageDataStatistics,vtkProcessObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Compute and return the average magnitude.
  double GetAverageMagnitude() {this->Update(); return this->AverageMagnitude;}

  // Description:
  // Compute and return the standard deviation.
  double GetStandardDeviation() {this->Update(); return this->StandardDeviation;}

  // Description:
  // Compute and return the standard deviation.
  long int GetCount() {this->Update(); return this->Count;}

  void Update();

  void SetInput(vtkImageData *input);
  vtkImageData *GetInput();


protected:
  vtkImageDataStatistics();
  ~vtkImageDataStatistics();
  vtkImageDataStatistics(const vtkImageDataStatistics&) {};
  void operator=(const vtkImageDataStatistics&) {};

  void Execute();

  double AverageMagnitude;
  double StandardDeviation;
  long int Count;

  vtkTimeStamp ExecuteTime;

};

#endif



