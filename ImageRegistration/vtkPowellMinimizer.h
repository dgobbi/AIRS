/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPowellMinimizer.h

=========================================================================*/
// .NAME vtkPowellMinimizer - use Powell's method to minimize a function
// .SECTION Description
// vtkPowellMinimizer will modify a set of parameters in order to find
// the minimum of a specified function.  This method conducts a series
// of linear searches and attempts to construct a conjugate set of search
// directions as it goes.

#ifndef __vtkPowellMinimizer_h
#define __vtkPowellMinimizer_h

#include "vtkObject.h"

class VTK_EXPORT vtkPowellMinimizer : public vtkObject
{
public:
  static vtkPowellMinimizer *New();
  vtkTypeMacro(vtkPowellMinimizer,vtkObject);
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
  virtual int Iterate();

  // Description:
  // Get the function value resulting from the minimization.
  vtkSetMacro(FunctionValue,double);
  double GetFunctionValue() { return this->FunctionValue; };

  // Description:
  // Set the amoeba contraction ratio.  The default value of 0.5 gives
  // fast convergence, but larger values such as 0.6 or 0.7 provide
  // greater stability.
  vtkSetClampMacro(ContractionRatio,double,0.5,1.0);
  vtkGetMacro(ContractionRatio,double);

  // Description:
  // Set the amoeba expansion ratio.  The default value is 2.0, which
  // provides rapid expansion.  Values between 1.1 and 2.0 are valid.
  vtkSetClampMacro(ExpansionRatio,double,1.0,2.0);
  vtkGetMacro(ExpansionRatio,double);

  // Description:
  // Specify the value tolerance to aim for during the minimization.
  vtkSetMacro(Tolerance,double);
  vtkGetMacro(Tolerance,double);

  // Description:
  // Specify the parameter tolerance to aim for during the minimization.
  vtkSetMacro(ParameterTolerance,double);
  vtkGetMacro(ParameterTolerance,double);

  // Description:
  // Specify the maximum number of iterations to try before giving up.
  vtkSetMacro(MaxIterations,int);
  vtkGetMacro(MaxIterations,int);

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
  vtkPowellMinimizer();
  ~vtkPowellMinimizer();

  void (*Function)(void *);
  void (*FunctionArgDelete)(void *);
  void *FunctionArg;

  int NumberOfParameters;
  char **ParameterNames;
  double *ParameterValues;
  double *ParameterScales;
  double FunctionValue;

  double ContractionRatio;
  double ExpansionRatio;

  double Tolerance;
  double ParameterTolerance;
  int MaxIterations;
  int Iterations;
  int FunctionEvaluations;

private:
  // Description:
  // Perform a golden-section search for a minimum.
  double PowellGolden(
    const double *p0, double y0, const double *v, double *p, int n,
    double gtol, bool *failed);

  // Description:
  // Initialize the workspace required for the method.
  void PowellInitialize();

  // Description:
  // Run one iteration of Powell's method.
  int PowellIterate();

  // Description:
  // Evaluate the function once.
  double PowellEvaluate(const double *p);

  double *PowellWorkspace;
  double **PowellVectors;

  vtkPowellMinimizer(const vtkPowellMinimizer&);  // Not implemented.
  void operator=(const vtkPowellMinimizer&);  // Not implemented.
};

#endif
