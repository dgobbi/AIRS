/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImageSingleMutualInformation.h,v $
  Language:  C++
  Date:      $Date: 2005/07/20 16:01:41 $
  Version:   $Revision: 1.5 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageSingleMutualInformation - 1-D Histogram for computing
// a Mutual Information histogram on two input images.
// .SECTION Description
// vtkImageSingleMutualInformation - This filter divides component space 
// into discrete bins.  It then counts the number of pixels associated
// with each bin.  The output is this "scatter plot".
// The dimensionality of the output depends on how many bins are wanted in 
// the image histograms.
// This filter can only handle images with 1 scalar component.
// The input can be any type. 
// The output is always int.
// Statistics are computed on the pixel values at the same time, and
// processed to provide a normalized mutual information value, in addition
// to image specific and inter-image entropy values.
// The SetStencilFunction, SetClippingExtents and ReverseStencil
// functions allow the statistics to be computed on an arbitrary
// portion of the input data.
// See the documentation for vtkImageStencil for more information.


#ifndef __vtkImageSingleMutualInformation_h
#define __vtkImageSingleMutualInformation_h

// #ifndef vtkFloatingPointType
// #define vtkFloatingPointType vtkFloatingPointType
// typedef float vtkFloatingPointType;
// #endif

#include "vtkImageToImageFilter.h"

class vtkImageStencilData;

class VTK_EXPORT vtkImageSingleMutualInformation : public vtkImageToImageFilter
{
public:
  static vtkImageSingleMutualInformation *New();
  vtkTypeRevisionMacro(vtkImageSingleMutualInformation,vtkImageToImageFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set/Get - The component spacing is the dimension of each bin.
  // This ends up being the spacing of the output "image".
  vtkSetMacro(ImageAComponentSpacing, float);
  vtkGetMacro(ImageAComponentSpacing, float);

  // Description:
  // Set/Get - The component origin is the location of bin (0, 0).
  // Note that if the Component extent does not include the value (0,0),
  // then this origin bin will not actually be in the output.
  // The origin of the output ends up being the same as the componenent origin.
  // For bins spanning the values 1000 to 2000,
  // this origin should be set to 1000.
  vtkSetMacro(ImageAComponentOrigin, float);
  vtkGetMacro(ImageAComponentOrigin, float);

  // Description:
  // Set/Get - The component extent sets the number/extent of the bins.
  // For an image histogram with 256 bins, the extent should be 0, 255.
  // The extent specifies inclusive min/max values.  
  // This implies the the top extent should be set to the number of bins - 1.
  void SetImageAComponentExtent(int extent[2]);
  void SetImageAComponentExtent(int min, int max);
  void GetImageAComponentExtent(int extent[2]);
  int *GetImageAComponentExtent() {return this->ImageAComponentExtent;}

  // Description:
  // Use a stencil to specify which voxels to include in the computation.
  void SetStencil(vtkImageStencilData *stencil);
  vtkImageStencilData *GetStencil();

  // Description:
  // Reverse the stencil.
  vtkSetMacro(ReverseStencil, int);
  vtkBooleanMacro(ReverseStencil, int);
  vtkGetMacro(ReverseStencil, int);

  // Description:
  // Get the normalized MI for the data.
  vtkGetMacro(NormalizedMI, double);

  // Description:
  // Get the mean value for the data.
  vtkGetMacro(MeanVoxel, double);
 
  
protected:
  vtkImageSingleMutualInformation();
  ~vtkImageSingleMutualInformation();

  float ImageAComponentSpacing;
  float ImageAComponentOrigin;
  int ImageAComponentExtent[2];

  void ExecuteInformation(vtkImageData **inputs, vtkImageData *output);
  void ComputeInputUpdateExtent(int inExt[6], int outExt[6], 
				int vtkNotUsed(whichInput));
  void ExecuteInformation(){this->vtkImageToImageFilter::ExecuteInformation();};
  void ExecuteData(vtkDataObject *out);

  double NormalizedMI;
  double MeanVoxel;

  int ReverseStencil;

private:
  vtkImageSingleMutualInformation(const vtkImageSingleMutualInformation&);
  // Not implemented.
  void operator=(const vtkImageSingleMutualInformation&);
  // Not implemented.
};

#endif



