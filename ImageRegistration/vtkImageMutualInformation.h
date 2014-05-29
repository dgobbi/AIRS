/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageMutualInformation.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageMutualInformation - Mutual information between two images
// .SECTION Description
// vtkImageMutualInformation generates a joint histogram from the two input
// images and uses this joint histogram to compute the mutual information
// between the images.  Each input image must have just a single scalar
// component.  The output of the filter will be the joint histogram, with
// the bins from the first input along the X axis, and the bins from the
// second image along the Y axis.  The number of bins and the bin size
// along each axis must be set before the filter executes.  After the
// filter has executed, the mutual information, normalized mutual information,
// and the value to minimize to register the images can be retrieved.
//
// References:
//
//  [1] D. Mattes, D.R. Haynor, H. Vesselle, T. Lewellen and W. Eubank,
//      PET-CT Image Registration in the Chest Using Free-form Deformations,
//      IEEE Transactions in Medical Imaging 22:120-128, 2003.
//
//  [2] C. Studholme, D.L.G. Hill and D.J. Hawkes,
//      An Overlap Invariant Measure of 3D Medical Image Alignment,
//      Pattern Recognition 32:71-86, 1999.

#ifndef __vtkImageMutualInformation_h
#define __vtkImageMutualInformation_h

#include "vtkThreadedImageAlgorithm.h"

class vtkImageStencilData;

class VTK_EXPORT vtkImageMutualInformation : public vtkThreadedImageAlgorithm
{
public:
  static vtkImageMutualInformation *New();
  vtkTypeMacro(vtkImageMutualInformation,vtkThreadedImageAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the type for the output.  The joint histogram will always be
  // computed using vtkIdType, but since vtkIdType is not directly
  // supported as an image data type, it will be converted to the requested
  // type for use as the output of the filter.  The default type is float.
  vtkSetMacro(OutputScalarType, int);
  vtkGetMacro(OutputScalarType, int);
  void SetOutputScalarTypeToFloat() {
    this->SetOutputScalarType(VTK_FLOAT); }
  void SetOutputScalarTypeToDouble() {
    this->SetOutputScalarType(VTK_DOUBLE); }
  void SetOutputScalarTypeToInt() {
    this->SetOutputScalarType(VTK_INT); }
  void SetOutputScalarTypeToUnsignedInt() {
    this->SetOutputScalarType(VTK_UNSIGNED_INT); }
  void SetOutputScalarTypeToLong() {
    this->SetOutputScalarType(VTK_LONG); }
  void SetOutputScalarTypeToUnsignedLong() {
    this->SetOutputScalarType(VTK_UNSIGNED_LONG); }
  void SetOutputScalarTypeToShort() {
    this->SetOutputScalarType(VTK_SHORT); }
  void SetOutputScalarTypeToUnsignedShort() {
    this->SetOutputScalarType(VTK_UNSIGNED_SHORT); }
  void SetOutputScalarTypeToSignedChar() {
    this->SetOutputScalarType(VTK_SIGNED_CHAR); }
  void SetOutputScalarTypeToUnsignedChar() {
    this->SetOutputScalarType(VTK_UNSIGNED_CHAR); }

  // Description:
  // Set the number of bins in the X and Y directions.
  vtkSetVector2Macro(NumberOfBins, int);
  vtkGetVector2Macro(NumberOfBins, int);

  // Description:
  // Set the center position of the first bin.  The default is zero.
  vtkSetVector2Macro(BinOrigin, double);
  vtkGetVector2Macro(BinOrigin, double);

  // Description:
  // Set the joint histogram bin spacing.  The default is one.
  vtkSetVector2Macro(BinSpacing, double);
  vtkGetVector2Macro(BinSpacing, double);

  // Description:
  // Use a stencil to limit the calculations to a specific region of
  // the input images.
  void SetStencilData(vtkImageStencilData *stencil);
  void SetStencil(vtkImageStencilData *stencil) {
    this->SetStencilData(stencil); }
  vtkImageStencilData *GetStencil();

  // Description:
  // Get the mutual information that was computed for the joint histogram.
  // The result is only valid after the filter has executed.
  vtkGetMacro(MutualInformation, double);

  // Description:
  // Get the normalized mutual information that was computed for the joint
  // histogram.  The result is only valid after the filter has executed.
  vtkGetMacro(NormalizedMutualInformation, double);

  // Description:
  // This is part of the executive, but is public so that it can be accessed
  // by non-member functions.
  virtual void ThreadedRequestData(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector,
                                   vtkImageData ***inData,
                                   vtkImageData **outData, int ext[6], int id);
protected:
  vtkImageMutualInformation();
  ~vtkImageMutualInformation();

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

  int NumberOfBins[2];
  double BinOrigin[2];
  double BinSpacing[2];

  int OutputScalarType;

  double MutualInformation;
  double NormalizedMutualInformation;

  vtkIdType *ThreadOutput[VTK_MAX_THREADS];
  int ThreadExecuted[VTK_MAX_THREADS];

private:
  vtkImageMutualInformation(const vtkImageMutualInformation&);  // Not implemented.
  void operator=(const vtkImageMutualInformation&);  // Not implemented.
};

#endif
