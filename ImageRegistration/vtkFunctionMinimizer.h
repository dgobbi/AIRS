/*=========================================================================

  Program:   AtamaiRegistration for VTK
  Module:    $RCSfile: vtkFunctionMinimizer.h,v $
  Language:  C++
  Date:      $Date: 2007/08/24 20:02:25 $
  Version:   $Revision: 1.4 $

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

// .NAME vtkFunctionMinimizer - nonlinear optimization with a simplex
// .SECTION Description
// vtkAmoebaMinimizer will modify a set of parameters in order to find
// the minimum of a specified function.  The method used is commonly
// known as the amoeba method, it constructs an n-dimensional simplex
// in parameter space (i.e. a tetrahedron if the number or parameters
// is 3) and moves the vertices around parameter space until a local
// minimum is found.  The amoeba method is robust, reasonably efficient,
// but is not guaranteed to find the global minimum if several local
// minima exist.
// .SECTION see also
// vtkAmoebaMinimizer vtkPowellMinimizer

#ifndef __vtkFunctionMinimizer_h
#define __vtkFunctionMinimizer_h

#include "vtkObject.h"

class VTK_EXPORT vtkFunctionMinimizer : public vtkObject
{
public:
  vtkTypeMacro(vtkFunctionMinimizer,vtkObject);
  static vtkFunctionMinimizer *New();

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify the function to be minimized.
  void SetFunction(void (*f)(void *), void *arg);

  // Description:
  // Set a function to call when a void* argument is being discarded.
  void SetFunctionArgDelete(void (*f)(void *));

  // Description:
  // Specify a variable to modify during the minimization.  Only the
  // variable you specify will be modified.  You must specify estimated
  // min and max possible values for each variable.  
  void SetScalarVariableBracket(const char *name, double min, double max);
  void SetScalarVariableBracket(const char *name, const double range[2]) {
    this->SetScalarVariableBracket(name,range[0],range[1]); };
  double *GetScalarVariableBracket(const char *name);
  void GetScalarVariableBracket(const char *name, double range[2]) {
    double *r = this->GetScalarVariableBracket(name);
    range[0] = r[0]; range[1] = r[1]; };

  // Description:
  // Get the value of a variable at the current stage of the minimization.
  double GetScalarVariableValue(const char *name);

  // Description:
  // Get a pointer to the scalar var array
  double *GetScalarVarPtr();

  // Description:
  // Iterate until the minimum is found to within the specified tolerance.
  void Minimize();

  // Description:
  // Initialize the minimization (this must be called before Iterate,
  // but is not necessary before Minimize).
  int Initialize();

  // Description:
  // Perform one iteration of minimization.
  // void Iterate();

  // Description:
  // Get the current value resulting from the minimization.
  // The SetScalarResult() method does not call Modified() on the
  // function minimizer (the Modified() call was causing performance
  // problems on multi-CPU machines)
  // vtkSetMacro(ScalarResult, double);
  void SetScalarResult(double result) { this->ScalarResult = result; }
  double GetScalarResult() { return this->ScalarResult; };

  // Description:
  // Specify the fractional tolerance to aim for during the minimization.
  vtkSetMacro(Tolerance,double);
  vtkGetMacro(Tolerance,double);

  // Description:
  // Specify the maximum number of iterations to try before 
  // printing an error and aborting.
  vtkSetMacro(MaxIterations,int);
  vtkGetMacro(MaxIterations,int);

  // Description:
  // Return the number of interations required for the last
  // minimization that was performed.
  vtkGetMacro(Iterations,int);

protected:
  vtkFunctionMinimizer();
  ~vtkFunctionMinimizer();

//BTX  
  void (*Function)(void *);
  void (*FunctionArgDelete)(void *);
  void *FunctionArg;
//ETX

  int NumberOfParameters;
  char **ParameterNames;
  int *ParameterIndices;
  double *Parameters;
//BTX
  double (*ParameterBrackets)[2];
//ETX
  double **Vertices;

  double ScalarResult;

  double Tolerance;
  int MaxIterations;
  int Iterations;

//BTX
  friend void vtkFunctionMinimizerFunction(void *data);
//ETX

private:
  vtkFunctionMinimizer(const vtkFunctionMinimizer&); // Not implemented.
  void operator=(const vtkFunctionMinimizer&); // Not implemented.

};

#endif
