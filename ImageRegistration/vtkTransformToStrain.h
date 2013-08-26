/*=========================================================================

  Program:   Atamai Image Registration and Segmentation
  Module:    vtkTransformToStrain.h

  Copyright (c) 2013 David Gobbi
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

  * Neither the name of David Gobbi, nor the names of any authors nor
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=========================================================================*/
// .NAME vtkTransformToStrain - create a strain map from a transform
// .SECTION Description
// vtkTransformToStrain takes any transform as input and produces a
// nine-component image where each voxel in the image is a strain
// tensor.  The Green's strain is used.

#ifndef __vtkTransformToStrain_h
#define __vtkTransformToStrain_h

#include "vtkAlgorithm.h"
#include "vtkImageData.h" // makes things a bit easier

class vtkAbstractTransform;

class VTK_HYBRID_EXPORT vtkTransformToStrain : public vtkAlgorithm
{
public:
  static vtkTransformToStrain *New();
  vtkTypeMacro(vtkTransformToStrain,vtkAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set/Get the transform which will be converted into a grid.
  virtual void SetInput(vtkAbstractTransform*);
  vtkGetObjectMacro(Input, vtkAbstractTransform);

  // Description:
  // Get/Set the extent of the grid.
  vtkSetVector6Macro(OutputExtent, int);
  vtkGetVector6Macro(OutputExtent, int);

  // Description:
  // Get/Set the origin of the grid.
  vtkSetVector3Macro(OutputOrigin, double);
  vtkGetVector3Macro(OutputOrigin, double);

  // Description:
  // Get/Set the spacing between samples in the grid.
  vtkSetVector3Macro(OutputSpacing, double);
  vtkGetVector3Macro(OutputSpacing, double);

  enum { GreensStrain, DeformationGradient };

  // Description:
  // Get/Set what to produce at the output.
  vtkSetMacro(OutputValue, int);
  vtkGetMacro(OutputValue, int);
  void SetOutputValueToGreensStrain() {
    this->SetOutputValue(GreensStrain); }
  void SetOutputValueToDeformationGradient() {
    this->SetOutputValue(DeformationGradient); }
  virtual const char *GetOutputValueAsString();

  // Description:
  // Get/Set the scalar type of the grid.  The default is float.
  vtkSetMacro(OutputScalarType, int);
  vtkGetMacro(OutputScalarType, int);
  void SetOutputScalarTypeToFloat() {
    this->SetOutputScalarType(VTK_DOUBLE); }
  void SetOutputScalarTypeToShort() {
    this->SetOutputScalarType(VTK_SHORT); }
  void SetOutputScalarTypeToUnsignedShort() {
    this->SetOutputScalarType(VTK_UNSIGNED_SHORT); }
  void SetOutputScalarTypeToUnsignedChar() {
    this->SetOutputScalarType(VTK_UNSIGNED_CHAR); }
  void SetOutputScalarTypeToSignedChar() {
    this->SetOutputScalarType(VTK_SIGNED_CHAR); }

  // Description:
  // Get the scale and shift to convert integer grid elements into
  // real values:  dx = scale*di + shift.  If the grid is of double type,
  // then scale = 1 and shift = 0.
  double GetValueScale() {
    this->UpdateShiftScale(); return this->ValueScale; };
  double GetValueShift() {
    this->UpdateShiftScale(); return this->ValueShift; };

  // Description:
  // Get the output data object for a port on this algorithm.
  vtkImageData* GetOutput();

  // Description:
  // see vtkAlgorithm for details
  virtual int ProcessRequest(
    vtkInformation*, vtkInformationVector**, vtkInformationVector*);

  // Description:
  // Compute the Green's strain from the deformation gradient tensor.
  static void ComputeGreensStrain(const double F[3][3], double G[3][3]);

  // Description:
  // Decompose a tensor into principal components, sortest from highest
  // to lowest.
  static void ComputePrincipals(
    const double F[3][3], double w[3], double G[3][3]);

protected:
  vtkTransformToStrain();
  ~vtkTransformToStrain();

  virtual void RequestInformation(
    vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  virtual void RequestData(
    vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  // Description:
  // Internal method to calculate the shift and scale values which
  // will provide maximum grid precision for a particular integer type.
  void UpdateShiftScale();

  unsigned long GetMTime();

  vtkAbstractTransform *Input;

  int OutputValue;
  int OutputScalarType;
  int OutputExtent[6];
  double OutputOrigin[3];
  double OutputSpacing[3];

  double ValueScale;
  double ValueShift;
  vtkTimeStamp ShiftScaleTime;

  // see algorithm for more info
  virtual int FillOutputPortInformation(int port, vtkInformation* info);

private:
  vtkTransformToStrain(const vtkTransformToStrain&);  // Not implemented.
  void operator=(const vtkTransformToStrain&);  // Not implemented.
};

#endif
