/*=========================================================================

  Module: vtkImageCorrelationRatio.h

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

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

#ifndef vtkImageCorrelationRatio_h
#define vtkImageCorrelationRatio_h

#include "vtkImageSimilarityMetric.h"

class vtkImageCorrelationRatioTLS;

class VTK_EXPORT vtkImageCorrelationRatio : public vtkImageSimilarityMetric
{
public:
  static vtkImageCorrelationRatio *New();
  vtkTypeMacro(vtkImageCorrelationRatio, vtkImageSimilarityMetric);

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the range of the data for the first input.
  // This is used to set the size of the array of partial sums that is used
  // to compute the metric.  The default range is (0, 255), which is only
  // suitable for 8-bit images.  For all other image types, a data range
  // must be provided.
  vtkSetVector2Macro(DataRange, double);

  // Description:
  // Get the correlation ratio between the two images.
  // The result is only valid after the filter has executed.
  vtkGetMacro(CorrelationRatio, double);

protected:
  vtkImageCorrelationRatio();
  ~vtkImageCorrelationRatio();

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

  double DataRange[2];

  int NumberOfBins;
  double BinOrigin;
  double BinSpacing;

  int OutputScalarType;

  double CorrelationRatio;

  vtkImageCorrelationRatioTLS *ThreadData;

private:
  vtkImageCorrelationRatio(const vtkImageCorrelationRatio&);  // Not implemented.
  void operator=(const vtkImageCorrelationRatio&);  // Not implemented.
};

#endif
