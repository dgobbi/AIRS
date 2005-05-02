/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImageFloodFill.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageFloodFill - Flood fill an image, given thresholds.
// .SECTION Description
// vtkImageFloodFill will perform a flood fill on an image, given upper
// and lower thresholds.  The flood filled region is the "inside" and
// will be filled with the original image by default, or with the
// InValue if ReplaceIn is set.  The "outside" is also filled with the
// original image values by default, or with OutValue if ReplaceOut
// is set.  Make sure that you have at least one of ReplaceIn or
// ReplaceOut set, or else the output will be the same as the imput.

#ifndef __vtkImageFloodFill_h
#define __vtkImageFloodFill_h

#include "vtkImageToImageFilter.h"
#include "vtkImageData.h"

class vtkPoints;
class vtkImageStencilData;

class VTK_EXPORT vtkImageFloodFill : public vtkImageToImageFilter
{
public:
  static vtkImageFloodFill *New();
  vtkTypeRevisionMacro(vtkImageFloodFill, vtkImageToImageFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the seeds.
  void SetSeedPoints(vtkPoints *points);
  vtkGetObjectMacro(SeedPoints, vtkPoints);

  // Description:
  // Limit the flood to a region (e.g. a slice).
  vtkSetVector6Macro(FloodExtent, int);
  vtkGetVector6Macro(FloodExtent, int);

  // Description:
  // The values greater than or equal to the value match.
  void ThresholdByUpper(double thresh);
  
  // Description:
  // The values less than or equal to the value match.
  void ThresholdByLower(double thresh);
  
  // Description:
  // The values in a range (inclusive) match
  void ThresholdBetween(double lower, double upper);

  // Description:
  // Determines whether to replace the pixel in range with InValue
  vtkSetMacro(ReplaceIn, int);
  vtkGetMacro(ReplaceIn, int);
  vtkBooleanMacro(ReplaceIn, int);
  
  // Description:
  // Replace the in range pixels with this value.
  void SetInValue(double val);
  vtkGetMacro(InValue, double);
  
  // Description:
  // Determines whether to replace the pixel out of range with OutValue
  vtkSetMacro(ReplaceOut, int);
  vtkGetMacro(ReplaceOut, int);
  vtkBooleanMacro(ReplaceOut, int);

  // Description:
  // Replace the in range pixels with this value.
  void SetOutValue(double val);
  vtkGetMacro(OutValue, double);
  
  // Description:
  // Get the Upper and Lower thresholds.
  vtkGetMacro(UpperThreshold, double);
  vtkGetMacro(LowerThreshold, double);
  
  // Description:
  // Specify the stencil to use.  The stencil can be created
  // from a vtkImplicitFunction or a vtkPolyData.
  virtual void SetStencil(vtkImageStencilData *stencil);
  vtkImageStencilData *GetStencil();

  // Description:
  // Reverse the stencil.
  vtkSetMacro(ReverseStencil, int);
  vtkBooleanMacro(ReverseStencil, int);
  vtkGetMacro(ReverseStencil, int);

  // Description:
  // Set the component to fill for multi-component images (default: 0)
  vtkSetMacro(ActiveComponent,int);
  vtkGetMacro(ActiveComponent,int);

  vtkSetObjectMacro(ImageMask, vtkImageData);
  vtkGetObjectMacro(ImageMask, vtkImageData);

  // Description:
  // Override the MTime to account for the seed points.
  unsigned long GetMTime();

  // Description:
  // Get the number of connected voxels in the flood fill.
  vtkSetMacro(NumberOfInVoxels, int);
  vtkGetMacro(NumberOfInVoxels, int);

  // Description:
  // Get the bounds of the flood fill region
  vtkSetVector6Macro(FloodBounds, int);
  vtkGetVector6Macro(FloodBounds, int);

protected:
  vtkImageFloodFill();
  ~vtkImageFloodFill();

  double UpperThreshold;
  double LowerThreshold;
  int ReplaceIn;
  double InValue;
  int ReplaceOut;
  double OutValue;
  
  vtkPoints *SeedPoints;
  int FloodExtent[6];
  int FloodBounds[6];

  int NumberOfInVoxels;

  int ReverseStencil;
  int ActiveComponent;

  vtkImageData *ImageMask;

  void ExecuteInformation(vtkImageData *inData, vtkImageData *outData);
  void ExecuteInformation(){this->vtkImageToImageFilter::ExecuteInformation();};
  void ExecuteData(vtkDataObject *out);

  void ComputeInputUpdateExtent(int inExt[6], int outExt[6]);

private:
  vtkImageFloodFill(const vtkImageFloodFill&);  // Not implemented.
  void operator=(const vtkImageFloodFill&);  // Not implemented.
};

#endif













