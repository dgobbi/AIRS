/*=========================================================================

  Module: vtkImageRegistration.h

  Copyright (c) 2006 Atamai, Inc.
  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageRegistration - Perform linear image registration.
// .SECTION Description
// This class will find the transformation that registers the source
// image to the target image.

#ifndef vtkImageRegistration_h
#define vtkImageRegistration_h

#include "vtkImageRegistrationBase.h"

class vtkImageStencilToImage;
class vtkImageToImageStencil;
class vtkImageReslice;
class vtkImageShiftScale;
class vtkImageSpread;
class vtkImageBSplineCoefficients;
class vtkFunctionMinimizer;
class vtkImageSimilarityMetric;

struct vtkImageRegistrationInfo;

class VTK_EXPORT vtkImageRegistration : public vtkImageRegistrationBase
{
public:
  vtkTypeMacro(vtkImageRegistration, vtkImageRegistrationBase);
  static vtkImageRegistration *New();
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Initialize the transform.  This will also initialize the
  // NumberOfEvaluations to zero.  If a TransformInitializer is
  // set, then only the rotation part of this matrix will be used,
  // and the initial translation will be set from the initializer.
  void Initialize(vtkMatrix4x4 *matrix);

  // Description:
  // Iterate the registration.  Returns zero if the termination condition has
  // been reached.
  int Iterate();

protected:
  vtkImageRegistration();
  ~vtkImageRegistration();

  int ExecuteRegistration();

  vtkMatrix4x4                    *InitialTransformMatrix;
  vtkImageReslice                 *ImageReslice;
  vtkImageBSplineCoefficients     *ImageBSpline;
  vtkImageShiftScale              *SourceImageTypecast;
  vtkImageShiftScale              *TargetImageTypecast;

  vtkImageSpread                  *Spread;
  vtkImageReslice                 *MaskReslice;
  vtkImageToImageStencil          *MaskToStencil;
  vtkImageStencilToImage          *StencilToMask;

  vtkImageRegistrationInfo        *RegistrationInfo;

private:
  // Copy constructor and assigment operator are purposely not implemented
  vtkImageRegistration(const vtkImageRegistration&);
  void operator=(const vtkImageRegistration&);
};

#endif /* vtkImageRegistration_h */
