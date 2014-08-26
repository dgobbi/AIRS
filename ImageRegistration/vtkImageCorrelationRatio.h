/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageCorrelationRatio.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageCorrelationRatio - Correlation ratio similarity metric.
// .SECTION Description
// vtkImageCorrelationRatio computes the correlation ratio for one image
// with respect to a second image.  Unlike many other image similarity
// metrics, it is not symmetrical. It assumes that there is a functional
// dependence of the first image on the second image, but does not assume
// that the reverse is true.  It is an efficient and robust method for
// multi-modal image registration.  For more information, please read
// the reference.
//
// References:
//
//  [1] A. Roche, G. Malandain, X. Pennec and N. Ayache,
//      The Correlation Ratio as a New Similarity Measure for Multimodal
//      Image Registration, MICCAI '98, LNCS 1496:1115-1124, 1998.

#ifndef __vtkImageCorrelationRatio_h
#define __vtkImageCorrelationRatio_h

#include "vtkThreadedImageAlgorithm.h"

class vtkImageStencilData;

class VTK_EXPORT vtkImageCorrelationRatio : public vtkThreadedImageAlgorithm
{
public:
  static vtkImageCorrelationRatio *New();
  vtkTypeMacro(vtkImageCorrelationRatio,vtkThreadedImageAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the range of the data for the first input.
  // This is used to set the size of the array of partial sums that is used
  // to compute the metric.  The default range is (0, 255), which is only
  // suitable for 8-bit images.  For all other image types, a data range
  // must be provided.
  vtkSetVector2Macro(DataRange, double);

  // Description:
  // Use a stencil to limit the calculations to a specific region of
  // the input images.
  void SetStencilData(vtkImageStencilData *stencil);
  void SetStencil(vtkImageStencilData *stencil) {
    this->SetStencilData(stencil); }
  vtkImageStencilData *GetStencil();

  // Description:
  // Get the correlation ratio between the two images.
  // The result is only valid after the filter has executed.
  vtkGetMacro(CorrelationRatio, double);

  // Description:
  // This is part of the executive, but is public so that it can be accessed
  // by non-member functions.
  virtual void ThreadedRequestData(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector,
                                   vtkImageData ***inData,
                                   vtkImageData **outData, int ext[6], int id);
protected:
  vtkImageCorrelationRatio();
  ~vtkImageCorrelationRatio();

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

  double DataRange[2];

  int NumberOfBins;
  double BinOrigin;
  double BinSpacing;

  int OutputScalarType;

  double CorrelationRatio;

  double *ThreadOutput[VTK_MAX_THREADS];
  bool ThreadExecuted[VTK_MAX_THREADS];

private:
  vtkImageCorrelationRatio(const vtkImageCorrelationRatio&);  // Not implemented.
  void operator=(const vtkImageCorrelationRatio&);  // Not implemented.
};

#endif
