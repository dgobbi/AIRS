/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageHistogramStatistics.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageHistogramStatistics - Compute statistics for an image
// .SECTION Description
// vtkImageHistogramStatistics computes statistics such as mean, median, and
// standard deviation.  These statistics are computed from the histogram
// of the image, rather than from the image itself, because this is more
// efficient than computing the statistics while traversing the pixels.
// If the input image is of type float or double, then the precision of
// will depend on what the MaximumNumberOfBins has been set to.
// .SECTION Thanks
// Thanks to David Gobbi at the Seaman Family MR Centre and Dept. of Clinical
// Neurosciences, Foothills Medical Centre, Calgary, for providing this class.

#ifndef __vtkImageHistogramStatistics_h
#define __vtkImageHistogramStatistics_h

#include "vtkImageHistogram.h"

class vtkImageStencilData;
class vtkIdTypeArray;

class VTK_EXPORT vtkImageHistogramStatistics : public vtkImageHistogram
{
public:
  static vtkImageHistogramStatistics *New();
  vtkTypeMacro(vtkImageHistogramStatistics,vtkImageHistogram);

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get the minimum value present in the image.  This value is computed
  // when Update() is called.
  double GetMinimum() { return this->Minimum; }

  // Description:
  // Get the maximum value present in the image.  This value is computed
  // when Update() is called.
  double GetMaximum() { return this->Maximum; }

  // Description:
  // Get the mean value of the image.  This value is computed when Update()
  // is called.
  double GetMean() { return this->Mean; }

  // Description:
  // Get the median value.  This is computed when Update() is called.
  double GetMedian() { return this->Median; }

  // Description:
  // Get the standard deviation of the values in the image.  This is
  // computed when Update() is called.
  double GetStandardDeviation() { return this->StandardDeviation; }

protected:
  vtkImageHistogramStatistics();
  ~vtkImageHistogramStatistics();

  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *);

  double Minimum;
  double Maximum;
  double Mean;
  double StandardDeviation;
  double Median;

private:
  vtkImageHistogramStatistics(const vtkImageHistogramStatistics&);  // Not implemented.
  void operator=(const vtkImageHistogramStatistics&);  // Not implemented.
};

#endif
