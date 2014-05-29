/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageSquaredDifference.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageSquaredDifference - Squared difference between images
// .SECTION Description
// vtkImageSquaredDifference computes the average squared difference of
// pixel values between two images. The images must have the same origin
// and spacing.

#ifndef __vtkImageSquaredDifference_h
#define __vtkImageSquaredDifference_h

#include "vtkThreadedImageAlgorithm.h"

class vtkImageStencilData;

class VTK_EXPORT vtkImageSquaredDifference : public vtkThreadedImageAlgorithm
{
public:
  static vtkImageSquaredDifference *New();
  vtkTypeMacro(vtkImageSquaredDifference,vtkThreadedImageAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Use a stencil to limit the calculations to a specific region of
  // the input images.
  void SetStencilData(vtkImageStencilData *stencil);
  void SetStencil(vtkImageStencilData *stencil) {
    this->SetStencilData(stencil); }
  vtkImageStencilData *GetStencil();

  // Description:
  // Get the cross correlation of the two images, with no normalization.
  // The result is only valid after the filter has executed.
  vtkGetMacro(SquaredDifference, double);

  // Description:
  // This is part of the executive, but is public so that it can be accessed
  // by non-member functions.
  virtual void ThreadedRequestData(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector,
                                   vtkImageData ***inData,
                                   vtkImageData **outData, int ext[6], int id);
protected:
  vtkImageSquaredDifference();
  ~vtkImageSquaredDifference();

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

  double SquaredDifference;

  double ThreadOutput[VTK_MAX_THREADS][2];

private:
  vtkImageSquaredDifference(const vtkImageSquaredDifference&);  // Not implemented.
  void operator=(const vtkImageSquaredDifference&);  // Not implemented.
};

#endif
