/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageHistogramStatistics.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageHistogramStatistics.h"

#include "vtkObjectFactory.h"
#include "vtkIdTypeArray.h"

#include <math.h>

vtkStandardNewMacro(vtkImageHistogramStatistics);

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageHistogramStatistics::vtkImageHistogramStatistics()
{
  this->AutomaticBinning = true;

  this->Minimum = 0;
  this->Maximum = 0;
  this->Median = 0;
  this->Mean = 0;
  this->StandardDeviation = 0;

  this->SetNumberOfOutputPorts(0);
}

//----------------------------------------------------------------------------
vtkImageHistogramStatistics::~vtkImageHistogramStatistics()
{
}

//----------------------------------------------------------------------------
void vtkImageHistogramStatistics::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Minimum: " << this->Minimum << "\n";
  os << indent << "Maximum: " << this->Maximum << "\n";
  os << indent << "Median: " << this->Median << "\n";
  os << indent << "Mean: " << this->Mean << "\n";
  os << indent << "StandardDeviation: " << this->StandardDeviation << "\n";
}

//----------------------------------------------------------------------------
int vtkImageHistogramStatistics::RequestData(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  this->Superclass::RequestData(request, inputVector, outputVector);

  vtkIdType total = this->Total;
  vtkIdType sum = 0;
  vtkIdType midSum = total/2;
  int midVal = 0;
  int minVal = -1;
  int maxVal = 0;
  double mom1 = 0;
  double mom2 = 0;
  int nx = this->Histogram->GetNumberOfTuples();
  vtkIdType *histogram = this->Histogram->GetPointer(0);
  for (int ix = 0; ix < nx; ++ix)
    {
    vtkIdType c = histogram[ix];
    sum += c;
    double dc = static_cast<double>(c);
    mom1 += dc*ix;
    mom2 += dc*ix*ix;
    midVal = (sum > midSum ? midVal : ix);
    minVal = (sum > 0 ? minVal : ix);
    maxVal = (c == 0 ? maxVal : ix);
    }
  if (minVal < maxVal)
    {
    minVal++;
    }

  double binSpacing = this->BinSpacing;
  double binOrigin = this->BinOrigin;

  // do the basic statistics
  this->Minimum = minVal*binSpacing + binOrigin;
  this->Maximum = maxVal*binSpacing + binOrigin;
  this->Median = midVal*binSpacing + binOrigin;
  this->Mean = 0.0;
  this->StandardDeviation = 0.0;
  if (total > 0)
    {
    this->Mean = mom1/total*binSpacing + binOrigin;
    }
  if (total > 1)
    {
    double term2 = mom1*mom1/total;
    if ((mom2 - term2) > 1e-10*mom2)
      {
      // use the fast method to compute standard deviation
      this->StandardDeviation = sqrt((mom2 - term2)/(total - 1))*binSpacing;
      }
    else
      {
      // use more accurate method to avoid cancellation error
      double xmean = mom1/total;
      for (int ix = 0; ix < nx; ++ix)
        {
        double ixd = xmean - ix;
        mom2 += ixd*ixd*histogram[ix];
        }
      this->StandardDeviation = sqrt(mom2/(total - 1))*binSpacing;
      }
    }

  return 1;
}
