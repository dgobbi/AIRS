/*=========================================================================

  Module: vtkImageSimilarityMetric.h

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/*! \class vtkImageSimilarityMetric
 *  \brief Base class for image similarity metrics.
 *
 *  This class exists to support the common paradigm for image registration
 *  that adjusts transformation parameters in order to optimize the measured
 *  similarity between two images.
 *
 *  The inputs of the metric are the two images to be compared and a stencil
 *  to mask the voxels that are to be compared.  The output of the metric is
 *  the value to minimize.
 */

#ifndef vtkImageSimilarityMetric_h
#define vtkImageSimilarityMetric_h

#include "vtkThreadedImageAlgorithm.h"

class vtkImageStencilData;
class vtkImageSimilarityMetricThreadData;
class vtkImageSimilarityMetricSMPThreadLocal;

class VTK_EXPORT vtkImageSimilarityMetric : public vtkThreadedImageAlgorithm
{
public:
  vtkTypeMacro(vtkImageSimilarityMetric,vtkThreadedImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  //@{
  //! Use a stencil to limit the comparison to a specific image region.
  /*!
   *  The stencil must mask out any voxels that are out-of-bounds due to
   *  the tranformation and resampling of either input image.
   */
  void SetStencilData(vtkImageStencilData *stencil);
  void SetStencil(vtkImageStencilData *stencil) {
    this->SetStencilData(stencil); }
  vtkImageStencilData *GetStencil();
  //@}

  //@{
  //! Set the estimated range of the specified input.
  /*!
   *  This range is used in a different manner for different metrics.
   *  For Mutual Information and other histogram-based metrics, it is
   *  the range of intensities that are to be mapped to histogram.
   *  For Cross Correlation, it is used to scale the reported cost value
   *  to a useful range between 0.0 and 1.0.  See the documentation for
   *  each metric to see how this range is used.
   */
  void SetInputRange(int idx, const double range[2]);
  void GetInputRange(int idx, double range[2]);
  //@}

  //@{
  //! Get the metric value.
  /*!
   *  This is the primary output of the metric, but it is not always the
   *  value that is used for image registration.  For registration, use
   *  GetCost() which returns a value that can be minimized.
   */
  double GetValue() { return this->Value; }

  //! Get the metric value in a form suitable for minimization.
  /*!
   *  This is the primary output of the metric.  In order to be able to
   *  use general purpose optimization code with metric, this must return
   *  a value that is smallest when the images are most similar.
   */
  double GetCost() { return this->Cost; }
  //@}

protected:
  vtkImageSimilarityMetric();
  ~vtkImageSimilarityMetric();

  int FillInputPortInformation(int port, vtkInformation *info);

  int RequestUpdateExtent(vtkInformation *request,
                          vtkInformationVector **inInfo,
                          vtkInformationVector *outInfo);

  //@{
  //! Subclasses override this method to create the thread-local info.
  /*!
   *  The subclass RequestData() method should construct the thread-local
   *  objects, and then call Superclass::RequestData(), i.e. this method,
   *  to cause PieceRequestData() and ReduceRequestData() to be called.
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inInfo,
                  vtkInformationVector *outInfo);

  //! Subclasses override this method to operate on one piece of input.
  /*!
   *  The whole extent of the input will be broken into non-overlapping
   *  pieces, and PieceRequestData() will be called for each piece.  This
   *  method will be called from multiple threads, and should use the
   *  thread-local object constructed in RequestData().
   */
  virtual void PieceRequestData(vtkInformation *request,
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector,
                                const int extent[6], vtkIdType piece) = 0;

  //! Subclasses override this method to reduce after multithreading.
  /*!
   *  After the threads have finished executing PieceRequestData() for all
   *  of the pieces, ReduceRequestData() will be called, and it should
   *  reduce the outputs for the various threads into a single output.
   *  This method should call SetCost().
   */
  virtual void ReduceRequestData(vtkInformation *request,
                                 vtkInformationVector **inInfo,
                                 vtkInformationVector *outInfo) = 0;

  //! Subclasses call this to set the metric value.
  void SetValue(double x) { this->Value = x; }

  //! Subclasses call this to set the value to be minimized.
  void SetCost(double x) { this->Cost = x; }
  //@}

private:
  vtkImageSimilarityMetric(const vtkImageSimilarityMetric&);
  void operator=(const vtkImageSimilarityMetric&);

  double InputRange[2][2];

  double Value;
  double Cost;

  friend class vtkImageSimilarityMetricFunctor;
  friend struct vtkImageSimilarityMetricThreadStruct;
};

#endif /* vtkImageSimilarityMetric_h */
