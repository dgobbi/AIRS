/*=========================================================================

  Module: vtkImageNeighborhoodCorrelation.h

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

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

#ifndef vtkImageNeighborhoodCorrelation_h
#define vtkImageNeighborhoodCorrelation_h

#include "vtkImageRegistrationModule.h" // For export macro
#include "vtkImageSimilarityMetric.h"

class vtkImageNeighborhoodCorrelationTLS;

class VTKIMAGEREGISTRATION_EXPORT vtkImageNeighborhoodCorrelation :
  public vtkImageSimilarityMetric
{
public:
  static vtkImageNeighborhoodCorrelation *New();
  vtkTypeMacro(vtkImageNeighborhoodCorrelation, vtkImageSimilarityMetric);

  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Description:
  // Set the neighborhood radius.  The neighborhood is a box function.
  // The default radius is 7, which gives a box of size 15*15*15.
  vtkSetVector3Macro(NeighborhoodRadius, int);
  vtkGetVector3Macro(NeighborhoodRadius, int);

protected:
  vtkImageNeighborhoodCorrelation();
  ~vtkImageNeighborhoodCorrelation();

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  void PieceRequestData(vtkInformation *request,
                        vtkInformationVector **inputVector,
                        vtkInformationVector *outputVector,
                        const int pieceExtent[6], vtkIdType pieceId) override;

  void ReduceRequestData(vtkInformation *request,
                         vtkInformationVector **inInfo,
                         vtkInformationVector *outInfo) override;

  int NeighborhoodRadius[3];

  vtkImageNeighborhoodCorrelationTLS *ThreadData;

private:
  vtkImageNeighborhoodCorrelation(const vtkImageNeighborhoodCorrelation&) = delete;
  void operator=(const vtkImageNeighborhoodCorrelation&) = delete;
};

#endif
