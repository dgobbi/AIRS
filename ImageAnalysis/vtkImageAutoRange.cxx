/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageAutoRange.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageAutoRange.h"

#include "vtkObjectFactory.h"
#include "vtkIdTypeArray.h"

#include <math.h>

vtkStandardNewMacro(vtkImageAutoRange);

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageAutoRange::vtkImageAutoRange()
{
  this->AutomaticBinning = true;

  this->AutoRange[0] = 0.0;
  this->AutoRange[1] = 1.0;

  this->FullRange[0] = 0.0;
  this->FullRange[1] = 0.0;

  this->SetNumberOfOutputPorts(0);
}

//----------------------------------------------------------------------------
vtkImageAutoRange::~vtkImageAutoRange()
{
}

//----------------------------------------------------------------------------
void vtkImageAutoRange::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "AutoRange: " << this->AutoRange[0] << " "
     << this->AutoRange[1] << "\n";
  os << indent << "FullRange: " << this->FullRange[0] << " "
     << this->FullRange[1] << "\n";
  os << indent << "AutoWindow: " << this->GetAutoWindow() << "\n";
  os << indent << "AutoLevel: " << this->GetAutoLevel() << "\n";
}

//----------------------------------------------------------------------------
int vtkImageAutoRange::RequestData(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  this->Superclass::RequestData(request, inputVector, outputVector);

  vtkIdType total = this->Total;
  vtkIdType sum = 0;
  vtkIdType lowSum = static_cast<vtkIdType>(total*0.005);
  vtkIdType highSum = static_cast<vtkIdType>(total*0.995);
  int lowVal = 0;
  int highVal = 0;
  int nextVal = 0;
  int minVal = -1;
  int maxVal = 0;
  int nx = this->Histogram->GetNumberOfTuples();
  vtkIdType *histogram = this->Histogram->GetPointer(0);
  for (int ix = 0; ix < nx; ++ix)
    {
    vtkIdType c = histogram[ix];
    sum += c;
    lowVal = (sum > lowSum ? lowVal : ix);
    highVal = (sum > highSum ? highVal : ix);
    nextVal = (sum > histogram[minVal+1] + 1 ? nextVal : ix);
    minVal = (sum > 0 ? minVal : ix);
    maxVal = (c == 0 ? maxVal : ix);
    }
  if (minVal < maxVal)
    {
    minVal++;
    }
  if (nextVal < maxVal)
    {
    nextVal++;
    }

  // do the autorange: first expand range by 10% at each end
  int e = static_cast<int>(0.10*(highVal - lowVal));
  lowVal -= e;
  highVal += e;

  // if there is an appreciable gap between the min value and
  // the next non-zero histogram bin, use nextVal as minimum
  // instead of lowVal
  //int d = static_cast<int>(0.10*(maxVal - minVal));
  //if (nextVal - minVal > d)
  //  {
  //  lowVal = nextVal;
  //  }

  double binSpacing = this->BinSpacing;
  double binOrigin = this->BinOrigin;

  this->AutoRange[0] = lowVal*binSpacing + binOrigin;
  this->AutoRange[1] = highVal*binSpacing + binOrigin;
  this->FullRange[0] = minVal*binSpacing + binOrigin;
  this->FullRange[1] = maxVal*binSpacing + binOrigin;
  // clamp the auto range to the full data range
  if (this->AutoRange[0] < this->FullRange[0])
    {
    this->AutoRange[0] = this->FullRange[0];
    }
  if (this->AutoRange[1] > this->FullRange[1])
    {
    this->AutoRange[1] = this->FullRange[1];
    }

  return 1;
}
