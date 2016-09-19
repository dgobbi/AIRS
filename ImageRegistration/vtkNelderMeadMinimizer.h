/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkNelderMeadMinimizer.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkNelderMeadMinimizer - nonlinear optimization with a simplex
// .SECTION Description
// vtkNelderMeadMinimizer will modify a set of parameters in order to find
// the minimum of a specified function.  The method used is commonly
// known as the amoeba method, it constructs an n-dimensional simplex
// in parameter space (i.e. a tetrahedron if the number or parameters
// is 3) and moves the vertices around parameter space until a local
// minimum is found.  The amoeba method is robust, reasonably efficient,
// but is not guaranteed to find the global minimum if several local
// minima exist.

#ifndef vtkNelderMeadMinimizer_h
#define vtkNelderMeadMinimizer_h

#include "vtkFunctionMinimizer.h"

class VTK_EXPORT vtkNelderMeadMinimizer : public vtkFunctionMinimizer
{
public:
  static vtkNelderMeadMinimizer *New();
  vtkTypeMacro(vtkNelderMeadMinimizer,vtkFunctionMinimizer);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the simplex contraction ratio.  The default value of 0.5 gives
  // fast convergence, but larger values such as 0.6 or 0.7 provide
  // greater stability.
  vtkSetClampMacro(ContractionRatio,double,0.5,1.0);
  vtkGetMacro(ContractionRatio,double);

  // Description:
  // Set the simplex expansion ratio.  The default value is 2.0, which
  // provides rapid expansion.  Values between 1.1 and 2.0 are valid.
  vtkSetClampMacro(ExpansionRatio,double,1.0,2.0);
  vtkGetMacro(ExpansionRatio,double);

protected:
  vtkNelderMeadMinimizer();
  ~vtkNelderMeadMinimizer();

  void Start();
  int Step();

  void InitializeAmoeba();
  void GetAmoebaParameterValues();
  void TerminateAmoeba();
  double TryAmoeba(double sum[], int high, double fac);
  int PerformAmoeba();
  int CheckParameterTolerance();

  double ContractionRatio;
  double ExpansionRatio;

  double **AmoebaVertices;
  double *AmoebaValues;
  double *AmoebaSum;
  double AmoebaSize;
  double AmoebaHighValue;
  int AmoebaNStepsNoImprovement;

private:
  vtkNelderMeadMinimizer(const vtkNelderMeadMinimizer&);  // Not implemented.
  void operator=(const vtkNelderMeadMinimizer&);  // Not implemented.
};

#endif
