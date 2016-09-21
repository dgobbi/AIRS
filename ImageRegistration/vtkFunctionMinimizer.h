/*=========================================================================

  Module: vtkFunctionMinimizer.h

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkFunctionMinimizer - base class for VTK function minimizers
// .SECTION Description
// This is the base class for methods that mininimize a cost function.

#ifndef vtkFunctionMinimizer_h
#define vtkFunctionMinimizer_h

#include "vtkObject.h"

class VTK_EXPORT vtkFunctionMinimizer : public vtkObject
{
public:
  vtkTypeMacro(vtkFunctionMinimizer,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify the function to be minimized.  When this function
  // is called, it must get the parameter values by calling
  // GetParameterValue() for each parameter, and then must
  // call SetFunctionValue() to tell the minimizer what the result
  // of the function evaluation was.  The number of function
  // evaluations used for the minimization can be retrieved
  // using GetFunctionEvaluations().
  void SetFunction(void (*f)(void *), void *arg);

  // Description:
  // Set a function to call when a void* argument is being discarded.
  void SetFunctionArgDelete(void (*f)(void *));

  // Description:
  // Set the initial value for the specified parameter.  Calling
  // this function for any parameter will reset the Iterations
  // and the FunctionEvaluations counts to zero.  You must also
  // use SetParameterScale() to specify the step size by which the
  // parameter will be modified during the minimization.  It is
  // preferable to specify parameters by name, rather than by
  // number.
  void SetParameterValue(const char *name, double value);
  void SetParameterValue(int i, double value);

  // Description:
  // Set the scale to use when modifying a parameter, i.e. the
  // initial amount by which the parameter will be modified
  // during the search for the minimum.  It is preferable to
  // identify scalars by name rather than by number.
  void SetParameterScale(const char *name, double scale);
  double GetParameterScale(const char *name);
  void SetParameterScale(int i, double scale);
  double GetParameterScale(int i) { return this->ParameterScales[i]; };

  // Description:
  // Get the value of a parameter at the current stage of the minimization.
  // Call this method within the function that you are minimizing in order
  // to get the current parameter values.  It is preferable to specify
  // parameters by name rather than by index.
  double GetParameterValue(const char *name);
  double GetParameterValue(int i) { return this->ParameterValues[i]; };

  // Description:
  // For completeness, an unchecked method to get the name for particular
  // parameter (the result will be NULL if no name was set).
  const char *GetParameterName(int i) { return this->ParameterNames[i]; };

  // Description:
  // Get the number of parameters that have been set.
  int GetNumberOfParameters() { return this->NumberOfParameters; };

  // Description:
  // Initialize the minimizer.  This will reset the number of parameters to
  // zero so that the minimizer can be reused.
  void Initialize();

  // Description:
  // Iterate until the minimum is found to within the specified tolerance,
  // or until the MaxIterations has been reached.
  virtual void Minimize();

  // Description:
  // Perform one iteration of minimization.  Returns zero if the tolerance
  // stopping criterion has been met.
  int Iterate();

  // Description:
  // Get the function value resulting from the minimization.
  vtkSetMacro(FunctionValue,double);
  double GetFunctionValue() { return this->FunctionValue; };

  // Description:
  // Specify the value tolerance to aim for during the minimization.
  // The minimizer continues iterating until (delta v) <= (tol*v).
  vtkSetMacro(Tolerance,double);
  vtkGetMacro(Tolerance,double);

  // Description:
  // Specify the parameter tolerance to aim for during the minimization.
  // The minimizer continues iterating until (delta p) <= tol*(p scale).
  vtkSetMacro(ParameterTolerance,double);
  vtkGetMacro(ParameterTolerance,double);

  // Description:
  // Specify the maximum number of iterations to try before giving up.
  vtkSetMacro(MaxIterations,int);
  vtkGetMacro(MaxIterations,int);

  // Description:
  // Tell the minimizer to abort before the next function evaluation.
  // This method can be called by the function that the minimizer is
  // calling, in order to tell the minimizer to return from Iterate()
  // or Minimize() before the next funtion call.
  void SetAbortFlag(int f) { this->AbortFlag = f; }
  int GetAbortFlag() { return this->AbortFlag; }
  void AbortFlagOn() { this->SetAbortFlag(1); }
  void AbortFlagOff() { this->SetAbortFlag(0); }

  // Description:
  // Return the number of interations that have been performed.  This
  // is not necessarily the same as the number of function evaluations.
  vtkGetMacro(Iterations,int);

  // Description:
  // Return the number of times that the function has been evaluated.
  vtkGetMacro(FunctionEvaluations,int);

  // Description:
  // Evaluate the function.  This is usually called internally by the
  // minimization code, but it is provided here as a public method.
  void EvaluateFunction();

protected:
  vtkFunctionMinimizer();
  ~vtkFunctionMinimizer();

  // Description:
  // Ask subclass to begin a new minimization.
  virtual void Start() = 0;

  // Description:
  // Ask subclass to perform one minimization step.
  virtual int Step() = 0;

  void (*Function)(void *);
  void (*FunctionArgDelete)(void *);
  void *FunctionArg;

  int NumberOfParameters;
  char **ParameterNames;
  double *ParameterValues;
  double *ParameterScales;
  double FunctionValue;

  double Tolerance;
  double ParameterTolerance;
  int MaxIterations;
  int Iterations;
  int FunctionEvaluations;
  int AbortFlag;

private:
  vtkFunctionMinimizer(const vtkFunctionMinimizer&);  // Not implemented.
  void operator=(const vtkFunctionMinimizer&);  // Not implemented.
};

#endif
