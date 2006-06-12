/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCenteredTransform.cxx,v $
  Language:  C++
  Date:      $Date: 2006/06/12 10:04:19 $
  Version:   $Revision: 1.3 $

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
#include "vtkCenteredTransform.h"

#include "vtkDataArray.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkMatrix4x4.h"

vtkCxxRevisionMacro(vtkCenteredTransform, "$Revision: 1.3 $");
vtkStandardNewMacro(vtkCenteredTransform);

//----------------------------------------------------------------------------
vtkCenteredTransform::vtkCenteredTransform()
{
  this->Center[0] = 0.0;
  this->Center[1] = 0.0;
  this->Center[2] = 0.0;

  this->Translation[0] = 0.0;
  this->Translation[1] = 0.0;
  this->Translation[2] = 0.0;
 
  this->RotationAnglesYXZ[0] = 0.0;
  this->RotationAnglesYXZ[1] = 0.0;
  this->RotationAnglesYXZ[2] = 0.0;
  
  this->IsotropicScale = 1.0;
}

//----------------------------------------------------------------------------
vtkCenteredTransform::~vtkCenteredTransform()
{
}

//----------------------------------------------------------------------------
void vtkCenteredTransform::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "Translation: " << "( "
     << this->Translation[0] << ", "
     << this->Translation[1] << ", "
     << this->Translation[2] << ")\n";

  os << indent << "Rotation Angle: " << "( "
     << this->RotationAnglesYXZ[0] << ", "
     << this->RotationAnglesYXZ[1] << ", "
     << this->RotationAnglesYXZ[2] << ")\n";
 
 os << indent << "Center: " << "( "
     << this->Center[0] << ", "
     << this->Center[1] << ", "
     << this->Center[2] << ")\n ";

 os << indent << "Scale: " << "( "
    << this->IsotropicScale << ")\n";
}

//----------------------------------------------------------------------------
void vtkCenteredTransformMatrixFromAngles(const double rotation[3],
                                          double matrix[3][3])
{ 
  double rx[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  double ry[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  double rz[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  double tmp1[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

  rx[1][1] = cos(rotation[0] * vtkMath::DoubleDegreesToRadians());
  rx[2][2] = cos(rotation[0] * vtkMath::DoubleDegreesToRadians());
  rx[1][2] = -sin(rotation[0] * vtkMath::DoubleDegreesToRadians());
  rx[2][1] = sin(rotation[0] * vtkMath::DoubleDegreesToRadians());
  rx[0][0] = 1.0;

  ry[0][0] = cos(rotation[1] * vtkMath::DoubleDegreesToRadians());
  ry[0][2] = sin(rotation[1] * vtkMath::DoubleDegreesToRadians());
  ry[2][0] = -sin(rotation[1] * vtkMath::DoubleDegreesToRadians());
  ry[2][2] = cos(rotation[1] * vtkMath::DoubleDegreesToRadians());
  ry[1][1] = 1.0;
   
  rz[0][0] = cos(rotation[2] * vtkMath::DoubleDegreesToRadians());
  rz[0][1] = -sin(rotation[2] * vtkMath::DoubleDegreesToRadians());
  rz[1][0] = sin(rotation[2] * vtkMath::DoubleDegreesToRadians());
  rz[1][1] = cos(rotation[2] * vtkMath::DoubleDegreesToRadians());
  rz[2][2] = 1.0;

  vtkMath::Multiply3x3(rx, ry, tmp1);
  vtkMath::Multiply3x3(rz, tmp1, matrix);
}

//----------------------------------------------------------------------------
void vtkCenteredTransformAnglesFromMatrix(const double matrix[3][3],
                                          double rotation[3])
{ 
  // use the 2nd and 3rd rows of the matrix to compute the angles
  double x2 = matrix[2][0];
  double y2 = matrix[2][1];
  double z2 = matrix[2][2];

  double x3 = matrix[1][0];
  double y3 = matrix[1][1];
  double z3 = matrix[1][2];

  // first find the rotation about the y axis
  double d1 = sqrt(x2*x2 + z2*z2);

  double cosTheta = z2/d1;
  double sinTheta = x2/d1;

  double theta = atan2(sinTheta, cosTheta);

  // now find rotation about x axis
  double d = sqrt(x2*x2 + y2*y2 + z2*z2);

  double sinPhi = y2/d;
  double cosPhi = (x2*x2 + z2*z2)/(d1*d);

  double phi = atan2(sinPhi, cosPhi);

  // finally, find rotation about z axis
  double x3p = x3*cosTheta - z3*sinTheta;
  double y3p = - sinPhi*sinTheta*x3 + cosPhi*y3 - sinPhi*cosTheta*z3;
  double d2 = sqrt(x3p*x3p + y3p*y3p);

  double cosAlpha = y3p/d2;
  double sinAlpha = x3p/d2;

  double alpha = atan2(sinAlpha, cosAlpha);

  // write out the result
  rotation[0] = phi/vtkMath::DoubleDegreesToRadians();
  rotation[1] = -theta/vtkMath::DoubleDegreesToRadians();
  rotation[2] = alpha/vtkMath::DoubleDegreesToRadians();
}

//----------------------------------------------------------------------------
// Update the 4x4 matrix. 
 void vtkCenteredTransform::InternalUpdate()
{
  this->Matrix->Identity();

  // do the computations to find the matrix

  // set the matrix values   
  double sc[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  double tmp2[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  double tmpvec[3] = {0,0,0};
  double final3x3[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  double centertranslation[3]= {0,0,0};
   
  centertranslation[0] = -this->Center[0];
  centertranslation[1] = -this->Center[1];
  centertranslation[2] = -this->Center[2];
  
  sc[0][0] = this->IsotropicScale;
  sc[1][1] = this->IsotropicScale;
  sc[2][2] = this->IsotropicScale;

  // compute the rotation matrix
  vtkCenteredTransformMatrixFromAngles(this->RotationAnglesYXZ, tmp2);
  
  // multiply by the scale factor
  vtkMath::Multiply3x3(tmp2, sc, final3x3);
  vtkMath::Multiply3x3(final3x3, centertranslation, tmpvec);

  // copy the rotation elements the to vtkMatrix4x4
  this->Matrix->Element[0][0] = final3x3[0][0];
  this->Matrix->Element[0][1] = final3x3[0][1];
  this->Matrix->Element[0][2] = final3x3[0][2];
  this->Matrix->Element[1][0] = final3x3[1][0];
  this->Matrix->Element[1][1] = final3x3[1][1];
  this->Matrix->Element[1][2] = final3x3[1][2];
  this->Matrix->Element[2][0] = final3x3[2][0];
  this->Matrix->Element[2][1] = final3x3[2][1];
  this->Matrix->Element[2][2] = final3x3[2][2];
  
  // copy the translation to the vtkMatrix4x4
  this->Matrix->Element[0][3] = 
    tmpvec[0] + this->Center[0] + this->Translation[0];
  this->Matrix->Element[1][3] =
    tmpvec[1] + this->Center[1] + this->Translation[1];
  this->Matrix->Element[2][3] =
    tmpvec[2] + this->Center[2] + this->Translation[2];

  // set the bottom row of the matrix to (0,0,0,1)
  this->Matrix->Element[3][0] = 0.0;
  this->Matrix->Element[3][1] = 0.0;
  this->Matrix->Element[3][2] = 0.0;
  this->Matrix->Element[3][3] = 1.0;
  
  this->Matrix->Modified();
} 

//------------------------------------------------------------------------
void vtkCenteredTransform::InternalDeepCopy(vtkAbstractTransform *transform)
{
  // copy all parameters to the new transform
  vtkCenteredTransform *t = (vtkCenteredTransform *)transform;
  this->Translation[0] = t->Translation[0]; 
  this->Translation[1] = t->Translation[1];
  this->Translation[2] = t->Translation[2];
 
  this->RotationAnglesYXZ[0] = t->RotationAnglesYXZ[0];
  this->RotationAnglesYXZ[1] = t->RotationAnglesYXZ[1];
  this->RotationAnglesYXZ[2] = t->RotationAnglesYXZ[2];

  this->IsotropicScale = t->IsotropicScale;

  this->Center[0] = t->Center[0];
  this->Center[1] = t->Center[1];
  this->Center[2] = t->Center[2];

  this->Modified();
}

//----------------------------------------------------------------------------
void vtkCenteredTransform::Inverse()
{
  // change all the parameters so that the transform goes in
  // the opposite direction
  double matrix[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  
  //Set the Center 
  this->Center[0] = this->Center[0] + this->Translation[0];
  this->Center[1] = this->Center[1] + this->Translation[1];
  this->Center[2] = this->Center[2] + this->Translation[2];

  this->Translation[0] = -this->Translation[0];
  this->Translation[1] = -this->Translation[1];
  this->Translation[2] = -this->Translation[2];

  if (this->IsotropicScale > 0)
    {
    this->IsotropicScale = 1/this->IsotropicScale;
    }
  else
    {
    vtkErrorMacro("IsotropicScale == 0, cannot invert the transform");
    }

  // create original rotation matrix
  vtkCenteredTransformMatrixFromAngles(this->RotationAnglesYXZ, matrix);
  // invert it by transposing it
  vtkMath::Transpose3x3(matrix, matrix);
  // get the new rotation angles
  vtkCenteredTransformAnglesFromMatrix(matrix, this->RotationAnglesYXZ);
  
  this->Modified();
  
  this->InternalUpdate();
}

//----------------------------------------------------------------------------
vtkAbstractTransform *vtkCenteredTransform::MakeTransform()
{
  return vtkCenteredTransform::New();
}



