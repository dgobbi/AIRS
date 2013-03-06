/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageNeighborhoodCorrelation.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageNeighborhoodCorrelation - Neighborhood-wise cross correlation
// .SECTION Description
// vtkImageNeighborhoodCorrelation computes the normalized cross
// correlation of two images, but does the normalization over small
// neighborhoods to make the metric robust to variations in signal
// intenstity across the image.  The two images must have the same
// origin and spacing, and must also have the same data type.

#ifndef __vtkImageNeighborhoodCorrelation_h
#define __vtkImageNeighborhoodCorrelation_h

#include "vtkThreadedImageAlgorithm.h"

class vtkImageStencilData;

class VTK_EXPORT vtkImageNeighborhoodCorrelation :
  public vtkThreadedImageAlgorithm
{
public:
  static vtkImageNeighborhoodCorrelation *New();
  vtkTypeMacro(vtkImageNeighborhoodCorrelation,vtkThreadedImageAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Use a stencil to limit the metric to a region of the images.
  void SetStencil(vtkImageStencilData *stencil);
  vtkImageStencilData *GetStencil();

  // Description:
  // Set the neighborhood radius.  The neighborhood is a box function.
  // The default radius is 4, which gives a box of size 9*9*9.
  vtkSetVector3Macro(NeighborhoodRadius, int);
  vtkGetVector3Macro(NeighborhoodRadius, int);

  // Description:
  // Get the metric value.
  vtkGetMacro(ValueToMinimize, double);

  // Description:
  // This is part of the executive, but is public so that it can be accessed
  // by non-member functions.
  virtual void ThreadedRequestData(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector,
                                   vtkImageData ***inData,
                                   vtkImageData **outData, int ext[6], int id);
protected:
  vtkImageNeighborhoodCorrelation();
  ~vtkImageNeighborhoodCorrelation();

  virtual int RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inInfo,
                                 vtkInformationVector *vtkNotUsed(outInfo));
  virtual int RequestInformation(vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inInfo,
                                 vtkInformationVector *vtkNotUsed(outInfo));
  virtual int RequestData(vtkInformation *,
			  vtkInformationVector **,
			  vtkInformationVector *);

  virtual int FillInputPortInformation(int port, vtkInformation *info);
  virtual int FillOutputPortInformation(int port, vtkInformation *info);

  int NeighborhoodRadius[3];
  double ValueToMinimize;
  double ThreadOutput[VTK_MAX_THREADS];

private:
  vtkImageNeighborhoodCorrelation(const vtkImageNeighborhoodCorrelation&);  // Not implemented.
  void operator=(const vtkImageNeighborhoodCorrelation&);  // Not implemented.
};

#endif
