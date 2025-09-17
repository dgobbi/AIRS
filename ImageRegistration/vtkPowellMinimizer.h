/*=========================================================================

  Module: vtkPowellMinimizer.h

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPowellMinimizer - use Powell's method to minimize a function
// .SECTION Description
// vtkPowellMinimizer will modify a set of parameters in order to find
// the minimum of a specified function.  This method conducts a series
// of linear searches and attempts to construct a conjugate set of search
// directions as it goes.

#ifndef vtkPowellMinimizer_h
#define vtkPowellMinimizer_h

#include "vtkImageRegistrationModule.h" // For export macro
#include "vtkFunctionMinimizer.h"

class VTKIMAGEREGISTRATION_EXPORT vtkPowellMinimizer :
  public vtkFunctionMinimizer
{
public:
  static vtkPowellMinimizer *New();
  vtkTypeMacro(vtkPowellMinimizer,vtkFunctionMinimizer);
  void PrintSelf(ostream& os, vtkIndent indent) override;

protected:
  vtkPowellMinimizer();
  ~vtkPowellMinimizer();

  // Description:
  // Use Brent's method to search for a minimum.
  double PowellBrent(
    const double *p0, double y0, const double *v, double *p, int n,
    const double bracket[3], double gtol);

  // Description:
  // Bracket a minimum, given a starting point.  Return the function
  // value at the center point of the bracket.
  double PowellBracket(
    const double *p0, double y0, const double *v, double *p, int n,
    double bracket[3], bool *failed);

  // Description:
  // Initialize the workspace required for the method.
  void Start() override;

  // Description:
  // Run one iteration of Powell's method.
  int Step() override;

  double *PowellWorkspace;
  double **PowellVectors;

private:
  vtkPowellMinimizer(const vtkPowellMinimizer&) = delete;
  void operator=(const vtkPowellMinimizer&) = delete;
};

#endif
