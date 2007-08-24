/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImageMutualInformation.h,v $
  Language:  C++
  Date:      $Date: 2007/08/24 20:02:25 $
  Version:   $Revision: 1.10 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageMutualInformation - 2 Dimensional Histogram for computing
// a Mutual Information 2D-histogram on two input images.
// .SECTION Description
// vtkImageMutualInformation - This filter divides component space for each
// image into discrete bins.  It then counts the number of pixels associated
// with each bin.  The output is this "scatter plot".
// The dimensionality of the output depends on how many bins are wanted in 
// the image histograms.
// This filter can only handle images with 1 scalar component.
// The input can be any type (but both must be the same type). 
// The output is always int.
// Statistics are computed on the pixel values at the same time, and
// processed to provide a normalized mutual information value, in addition
// to image specific and inter-image entropy values.
// The SetStencilFunction, SetClippingExtents and ReverseStencil
// functions allow the statistics to be computed on an arbitrary
// portion of the input data.
// See the documentation for vtkImageStencil for more information.


#ifndef __vtkImageMutualInformation_h
#define __vtkImageMutualInformation_h

#include "vtkSystemIncludes.h"
#ifndef vtkFloatingPointType
#define vtkFloatingPointType vtkFloatingPointType
typedef float vtkFloatingPointType;
#endif

#include "vtkImageTwoInputFilter.h"

class vtkImageStencilData;
class vtkInformation;

class VTK_EXPORT vtkImageMutualInformation : public vtkImageTwoInputFilter
{
public:
  vtkTypeRevisionMacro(vtkImageMutualInformation,vtkImageTwoInputFilter);
  static vtkImageMutualInformation *New();

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set/Get - The component spacing is the dimension of each bin.
  // This ends up being the spacing of the output "image".
  vtkSetMacro(ImageAComponentSpacing, float);
  vtkGetMacro(ImageAComponentSpacing, float);
  vtkSetMacro(ImageBComponentSpacing, float);
  vtkGetMacro(ImageBComponentSpacing, float);

  // Description:
  // Set/Get - The component origin is the location of bin (0, 0).
  // Note that if the Component extent does not include the value (0,0),
  // then this origin bin will not actually be in the output.
  // The origin of the output ends up being the same as the componenent origin.
  // For bins spanning the values 1000 to 2000,
  // this origin should be set to 1000.
  vtkSetMacro(ImageAComponentOrigin, float);
  vtkGetMacro(ImageAComponentOrigin, float);
  vtkSetMacro(ImageBComponentOrigin, float);
  vtkGetMacro(ImageBComponentOrigin, float);

  // Description:
  // Set/Get - The component extent sets the number/extent of the bins.
  // For an image histogram with 256 bins, the extent should be 0, 255.
  // The extent specifies inclusive min/max values.  
  // This implies the the top extent should be set to the number of bins - 1.
  void SetImageAComponentExtent(int extent[2]);
  void SetImageAComponentExtent(int min, int max);
  void GetImageAComponentExtent(int extent[2]);
  int *GetImageAComponentExtent() {return this->ImageAComponentExtent;}
  void SetImageBComponentExtent(int extent[2]);
  void SetImageBComponentExtent(int min, int max);
  void GetImageBComponentExtent(int extent[2]);
  int *GetImageBComponentExtent() {return this->ImageBComponentExtent;}


  // Description:
  // Override superclass setinput1 methods.
  virtual void SetInput1(vtkImageData *input);
  virtual vtkImageData *GetInput1();

  // Description:
  // Override superclass setinput2 methods.
  virtual void SetInput2(vtkImageData *input);
  virtual vtkImageData *GetInput2();

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
 
  
protected:
  vtkImageMutualInformation();
  ~vtkImageMutualInformation();

#if (VTK_MAJOR_VERSION >= 5) 
  // see vtkAlgorithm for docs.
  virtual int FillInputPortInformation(int, vtkInformation*);
#endif

  void ExecuteInformation(vtkImageData **inputs, vtkImageData *output);
  void ComputeInputUpdateExtent(int inExt[6], int outExt[6], 
				int vtkNotUsed(whichInput));
  void ExecuteInformation(){this->vtkImageTwoInputFilter::ExecuteInformation();};
  void ExecuteData(vtkDataObject *out);

  float ImageAComponentSpacing;
  float ImageAComponentOrigin;
  int ImageAComponentExtent[2];
  float ImageBComponentSpacing;
  float ImageBComponentOrigin;
  int ImageBComponentExtent[2];

  double NormalizedMI;
  int ReverseStencil;

private:
  vtkImageMutualInformation(const vtkImageMutualInformation&);  // Not implemented.
  void operator=(const vtkImageMutualInformation&);  // Not implemented.
};

#endif



