/*=========================================================================

  Module: vtkImageCrossCorrelation.h

  Copyright (c) 2006 Atamai, Inc.
  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageCrossCorrelation - Cross correlation between two images
// .SECTION Description
// vtkImageCrossCorrelation computes the cross correlation and the normalized
// cross correlation of two input images.  The images must have the same
// origin and spacing.

#ifndef vtkImageCrossCorrelation_h
#define vtkImageCrossCorrelation_h

#include "vtkImageRegistrationModule.h" // For export macro
#include "vtkImageSimilarityMetric.h"

class vtkImageCrossCorrelationTLS;

class VTKIMAGEREGISTRATION_EXPORT vtkImageCrossCorrelation :
  public vtkImageSimilarityMetric
{
public:
  static vtkImageCrossCorrelation *New();
  vtkTypeMacro(vtkImageCrossCorrelation, vtkImageSimilarityMetric);

  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Description:
  // Get the cross correlation of the two images, with no normalization.
  // The result is only valid after the filter has executed.
  vtkGetMacro(CrossCorrelation, double);

  // Description:
  // Get the normalized cross correlation of the two images.
  // The result is only valid after the filter has executed.
  vtkGetMacro(NormalizedCrossCorrelation, double);

  // The metrics (CrossCorrelation, NormalizedCrossCorrelation).
  enum { CC, NCC };

  // Description:
  // Set the metric to use for the cost. The default metric is Normalized
  // Cross Correlation.
  void SetMetricToCrossCorrelation() { this->SetMetric(CC); }
  void SetMetricToNormalizedCrossCorrelation() { this->SetMetric(NCC); }
  vtkSetMacro(Metric, int);
  vtkGetMacro(Metric, int);

protected:
  vtkImageCrossCorrelation();
  ~vtkImageCrossCorrelation();

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

  int Metric;

  double CrossCorrelation;
  double NormalizedCrossCorrelation;

  vtkImageCrossCorrelationTLS *ThreadData;

private:
  vtkImageCrossCorrelation(const vtkImageCrossCorrelation&) = delete;
  void operator=(const vtkImageCrossCorrelation&) = delete;
};

#endif
