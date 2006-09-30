/*=========================================================================

  Program:   AtamaiRegistration for VTK
  Module:    $RCSfile: vtkCenteredTransform.h,v $
  Language:  C++
  Date:      $Date: 2006/09/30 02:46:34 $
  Version:   $Revision: 1.2 $

Copyright (c) 2006 Atamai, Inc.
All rights reserved.

THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE SOFTWARE "AS IS"
WITHOUT EXPRESSED OR IMPLIED WARRANTY INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE.  IN NO EVENT SHALL ANY COPYRIGHT HOLDER OR OTHER PARTY WHO MAY
MODIFY AND/OR REDISTRIBUTE THE SOFTWARE UNDER THE TERMS OF THIS LICENSE
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, LOSS OF DATA OR DATA BECOMING INACCURATE
OR LOSS OF PROFIT OR BUSINESS INTERRUPTION) ARISING IN ANY WAY OUT OF
THE USE OR INABILITY TO USE THE SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES.
=========================================================================*/

// .NAME vtkCenteredTransform - a transform that performs the
// transformation around Center instead of the Origin
// .SECTION Description
// vtkCenteredTransform is a transformation which will compute a 4X4
// Homogeneous transformation matrix for the operations Scale, Rotation and
// Translation respectively with respect to the user selected
// Center. The transformation performs the operation in the following
// order (a) Translate(-Center), (b) Scale(Scale) (c) RotateY(Ry)
// (d)RotateX(Rx) (e) RotateZ(Rz) (f) Translate(Center) (g)
// Translate(Translation). Note that the Scale is Isotropic. 
// .SECTION see also
// vtkLinearTransform


#ifndef __vtkCenteredTransform_h
#define __vtkCenteredTransform_h

#include "vtkLinearTransform.h"

class VTK_EXPORT vtkCenteredTransform : public vtkLinearTransform
{
public:
  static vtkCenteredTransform *New();

  vtkTypeRevisionMacro(vtkCenteredTransform,vtkLinearTransform);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the center of the centered transform.
  vtkSetVector3Macro(Center, double);
  vtkGetVector3Macro(Center, double);

  // Description:
  // Set the translation for the centered transform.
  vtkSetVector3Macro(Translation, double);
  vtkGetVector3Macro(Translation, double);

  // Description:
  // Set the rotation angles for the centered transform.
  vtkSetVector3Macro(RotationAnglesYXZ, double);
  vtkGetVector3Macro(RotationAnglesYXZ, double);

  // Description:
  // Set the scale for the centered transform.
  vtkSetMacro(IsotropicScale, double);
  vtkGetMacro(IsotropicScale, double);

  // Description:
  // Initialize the transform, given a matrix and
  // a value to use as the center.
  void Initialize(vtkMatrix4x4 *matrix, double center[3]);

  // Description:
  // Inverts the transformation.  This method will change the values
  // of all the parameters so that the transform becomes inverted.
  void Inverse();

  // Description:
  // Create another transform with all parameters set to default values.
  vtkAbstractTransform *MakeTransform();

protected:
  vtkCenteredTransform();
  ~vtkCenteredTransform();

  // Description:
  // Update the matrix from the values of the parameters.
  void InternalUpdate();

  // Description:
  // Copy all parameters from the other transform.
  void InternalDeepCopy(vtkAbstractTransform *t);

    
  double Translation[3];
  double RotationAnglesYXZ[3];
  double IsotropicScale;
  double Center[3];

private:
  vtkCenteredTransform(const vtkCenteredTransform&);  // Not implemented.
  void operator=(const vtkCenteredTransform&);  // Not implemented.
};

#endif





