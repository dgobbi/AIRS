/*=========================================================================

  Module: vtkImageSquaredDifference.h

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageSquaredDifference - Squared difference between images
// .SECTION Description
// vtkImageSquaredDifference computes the average squared difference of
// pixel values between two images. The images must have the same origin
// and spacing.

#ifndef vtkImageSquaredDifference_h
#define vtkImageSquaredDifference_h

#include "vtkImageSimilarityMetric.h"

class vtkImageSquaredDifferenceTLS;

class VTK_EXPORT vtkImageSquaredDifference : public vtkImageSimilarityMetric
{
public:
  static vtkImageSquaredDifference *New();
  vtkTypeMacro(vtkImageSquaredDifference, vtkImageSimilarityMetric);

  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkImageSquaredDifference();
  ~vtkImageSquaredDifference();

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector);

  void PieceRequestData(vtkInformation *request,
                        vtkInformationVector **inputVector,
                        vtkInformationVector *outputVector,
                        const int pieceExtent[6], vtkIdType pieceId);

  void ReduceRequestData(vtkInformation *request,
                         vtkInformationVector **inInfo,
                         vtkInformationVector *outInfo);

  vtkImageSquaredDifferenceTLS *ThreadData;

private:
  vtkImageSquaredDifference(const vtkImageSquaredDifference&);  // Not implemented.
  void operator=(const vtkImageSquaredDifference&);  // Not implemented.
};

#endif
