/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageStatistics.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageStatistics - Compute statistics for an image
// .SECTION Description
// vtkImageStatistics computes statistics such as mean, median, and
// standard deviation.  These statistics are computed from the histogram
// of the image, rather than from the image itself, because this is more
// efficient for large images.

#ifndef __vtkImageStatistics_h
#define __vtkImageStatistics_h

#include "vtkImageHistogram.h"

class vtkImageStencilData;
class vtkIdTypeArray;

class VTK_EXPORT vtkImageStatistics : public vtkImageHistogram
{
public:
  static vtkImageStatistics *New();
  vtkTypeMacro(vtkImageStatistics,vtkImageHistogram);

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get basic statistics, these are computed from the histogram instead
  // of being computed directly from the data in order to make the
  // computation more efficient for large image.  Update must be called
  // before this information is retrieved.
  double GetMinimum() { return this->Minimum; }
  double GetMaximum() { return this->Maximum; }
  double GetMean() { return this->Mean; }
  double GetMedian() { return this->Median; }
  double GetStandardDeviation() { return this->StandardDeviation; }

protected:
  vtkImageStatistics();
  ~vtkImageStatistics();

  virtual void ComputeStatistics();

  double Minimum;
  double Maximum;
  double Mean;
  double StandardDeviation;
  double Median;

private:
  vtkImageStatistics(const vtkImageStatistics&);  // Not implemented.
  void operator=(const vtkImageStatistics&);  // Not implemented.
};

#endif
