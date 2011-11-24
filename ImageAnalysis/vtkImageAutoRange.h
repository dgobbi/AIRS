/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageAutoRange.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageAutoRange - Compute a range for viewing an image
// .SECTION Description
// vtkImageAutoRange computes window/level parameters for viewing an
// image by ignoring voxels that are much brighter than the rest of the
// image.  This provides better results than using the full range of
// the data, which often results in the image looking too dark.
// .SECTION Thanks
// Thanks to David Gobbi at the Seaman Family MR Centre and Dept. of Clinical
// Neurosciences, Foothills Medical Centre, Calgary, for providing this class.

#ifndef __vtkImageAutoRange_h
#define __vtkImageAutoRange_h

#include "vtkImageHistogram.h"

class vtkImageStencilData;
class vtkIdTypeArray;

class VTK_EXPORT vtkImageAutoRange : public vtkImageHistogram
{
public:
  static vtkImageAutoRange *New();
  vtkTypeMacro(vtkImageAutoRange,vtkImageHistogram);

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get an automatically computed range for use by lookup tables.  This
  // range will exclude outliers that are beyond the 99th percentile, meaning
  // that the brightest and darkest two percent of the voxels will be excluded
  // unless they are only 10 percent brighter (or darker) than the 99th
  // percentile.  This helps to avoid showing the images as overly dark just
  // because a few voxels are extremely bright compared to the rest.  You
  // must call Update before calling this method.
  vtkGetVector2Macro(AutoRange, double);

  // Description:
  // Get automatically computed ColorWindow and ColorLevel values for use
  // by image viewers.  These are computed in the same manner as AutoRange.
  // You must call Update before calling these methods.
  double GetAutoWindow() { return (this->AutoRange[1]-this->AutoRange[0]); }
  double GetAutoLevel() { return 0.5*(this->AutoRange[0]+this->AutoRange[1]); }

  // Description:
  // Get the full range of the data.  You must call Update before calling
  // this method.
  vtkGetVector2Macro(FullRange, double);

protected:
  vtkImageAutoRange();
  ~vtkImageAutoRange();

  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *);

  double AutoRange[2];
  double FullRange[2];

private:
  vtkImageAutoRange(const vtkImageAutoRange&);  // Not implemented.
  void operator=(const vtkImageAutoRange&);  // Not implemented.
};

#endif
