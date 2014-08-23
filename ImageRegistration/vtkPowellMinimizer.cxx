/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPowellMinimizer.cxx

=========================================================================*/
#include "vtkPowellMinimizer.h"
#include "vtkObjectFactory.h"

#define  N_STEPS_NO_VALUE_IMPROVEMENT  2
#define  N_STEPS_NO_PARAM_IMPROVEMENT  18

vtkStandardNewMacro(vtkPowellMinimizer);

//----------------------------------------------------------------------------
vtkPowellMinimizer::vtkPowellMinimizer()
{
  this->Function = NULL;
  this->FunctionArg = NULL;
  this->FunctionArgDelete = NULL;

  this->NumberOfParameters = 0;
  this->ParameterNames = NULL;
  this->ParameterValues = NULL;
  this->ParameterScales = NULL;

  this->FunctionValue = 0.0;

  this->ContractionRatio = 0.5;
  this->ExpansionRatio = 2.0;

  this->Tolerance = 1e-4;
  this->ParameterTolerance = 1e-4;
  this->MaxIterations = 1000;
  this->Iterations = 0;
  this->FunctionEvaluations = 0;

  // specific to Powell's method
  this->PowellWorkspace = 0;
  this->PowellVectors = 0;
}

//----------------------------------------------------------------------------
vtkPowellMinimizer::~vtkPowellMinimizer()
{
  if ((this->FunctionArg) && (this->FunctionArgDelete))
    {
    (*this->FunctionArgDelete)(this->FunctionArg);
    }
  this->FunctionArg = NULL;
  this->FunctionArgDelete = NULL;
  this->Function = NULL;

  if (this->ParameterNames)
    {
    for (int i = 0; i < this->NumberOfParameters; i++)
      {
      if (this->ParameterNames[i])
        {
        delete [] this->ParameterNames[i];
        }
      }
    delete [] this->ParameterNames;
    this->ParameterNames = NULL;
    }
  if (this->ParameterValues)
    {
    delete [] this->ParameterValues;
    this->ParameterValues = NULL;
    }
  if (this->ParameterScales)
    {
    delete [] this->ParameterScales;
    this->ParameterScales = NULL;
    }

  this->NumberOfParameters = 0;

  // specific to Powell's method
  delete [] this->PowellVectors;
  delete [] this->PowellWorkspace;
}

//----------------------------------------------------------------------------
void vtkPowellMinimizer::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "NumberOfParameters: " << this->GetNumberOfParameters() << "\n";
  if (this->NumberOfParameters > 0)
    {
    int i;

    os << indent << "ParameterValues: \n";
    for (i = 0; i < this->NumberOfParameters; i++)
      {
      const char *name = this->GetParameterName(i);
      os << indent << "  ";
      if (name)
        {
        os << name << ": ";
        }
      else
        {
        os << i << ": ";
        }
      os << this->GetParameterValue(i) << "\n";
      }

    os << indent << "ParameterScales: \n";
    for (i = 0; i < this->NumberOfParameters; i++)
      {
      const char *name = this->GetParameterName(i);
      os << indent << "  ";
      if (name)
        {
        os << name << ": ";
        }
      else
        {
        os << i << ": ";
        }
      os << this->GetParameterScale(i) << "\n";
      }
    }

  os << indent << "FunctionValue: " << this->GetFunctionValue() << "\n";
  os << indent << "FunctionEvaluations: " << this->GetFunctionEvaluations()
     << "\n";
  os << indent << "Iterations: " << this->GetIterations() << "\n";
  os << indent << "MaxIterations: " << this->GetMaxIterations() << "\n";
  os << indent << "Tolerance: " << this->GetTolerance() << "\n";
  os << indent << "ParameterTolerance: " << this->GetParameterTolerance() << "\n";
  os << indent << "ContractionRatio: " << this->GetContractionRatio() << "\n";
  os << indent << "ExpansionRatio: " << this->GetExpansionRatio() << "\n";
}

//----------------------------------------------------------------------------
void vtkPowellMinimizer::SetFunction(void (*f)(void *), void *arg)
{
  if ( f != this->Function || arg != this->FunctionArg )
    {
    // delete the current arg if there is one and a delete meth
    if ((this->FunctionArg) && (this->FunctionArgDelete))
      {
      (*this->FunctionArgDelete)(this->FunctionArg);
      }
    this->Function = f;
    this->FunctionArg = arg;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkPowellMinimizer::SetFunctionArgDelete(void (*f)(void *))
{
  if ( f != this->FunctionArgDelete)
    {
    this->FunctionArgDelete = f;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
double vtkPowellMinimizer::GetParameterValue(const char *name)
{
  for (int i = 0; i < this->NumberOfParameters; i++)
    {
    if (this->ParameterNames[i] && strcmp(name,this->ParameterNames[i]) == 0)
      {
      return this->ParameterValues[i];
      }
    }
  vtkErrorMacro("GetParameterValue: no parameter named " << name);
  return 0.0;
}

//----------------------------------------------------------------------------
void vtkPowellMinimizer::SetParameterValue(const char *name, double val)
{
  int i;

  for (i = 0; i < this->NumberOfParameters; i++)
    {
    if (this->ParameterNames[i] && strcmp(name,this->ParameterNames[i]) == 0)
      {
      break;
      }
    }

  this->SetParameterValue(i, val);

  if (!this->ParameterNames[i])
    {
    char *cp = new char[strlen(name)+8];
    strcpy(cp,name);
    this->ParameterNames[i] = cp;
    }
}

//----------------------------------------------------------------------------
void vtkPowellMinimizer::SetParameterValue(int i, double val)
{
  if (i < this->NumberOfParameters)
    {
    if (this->ParameterValues[i] != val)
      {
      this->ParameterValues[i] = val;
      this->Iterations = 0; // reset to start
      this->FunctionEvaluations = 0;
      this->Modified();
      }
    return;
    }

  int n = this->NumberOfParameters + 1;

  char **newParameterNames = new char *[n];
  double *newParameterValues = new double[n];
  double *newParameterScales = new double[n];

  for (int j = 0; j < this->NumberOfParameters; j++)
    {
    newParameterNames[j] = this->ParameterNames[j];
    this->ParameterNames[j] = NULL; // or else it will be deleted in Initialize
    newParameterValues[j] = this->ParameterValues[j];
    newParameterScales[j] = this->ParameterScales[j];
    }

  newParameterNames[n-1] = 0;
  newParameterValues[n-1] = val;
  newParameterScales[n-1] = 1.0;

  this->Initialize();

  this->NumberOfParameters = n;
  this->ParameterNames = newParameterNames;
  this->ParameterValues = newParameterValues;
  this->ParameterScales = newParameterScales;

  this->Iterations = 0; // reset to start
  this->FunctionEvaluations = 0;
}

//----------------------------------------------------------------------------
double vtkPowellMinimizer::GetParameterScale(const char *name)
{
  for (int i = 0; i < this->NumberOfParameters; i++)
    {
    if (this->ParameterNames[i] && strcmp(name,this->ParameterNames[i]) == 0)
      {
      return this->ParameterScales[i];
      }
    }
  vtkErrorMacro("GetParameterScale: no parameter named " << name);
  return 1.0;
}

//----------------------------------------------------------------------------
void vtkPowellMinimizer::SetParameterScale(const char *name, double scale)
{
  for (int i = 0; i < this->NumberOfParameters; i++)
    {
    if (this->ParameterNames[i] && strcmp(name,this->ParameterNames[i]) == 0)
      {
      this->SetParameterScale(i, scale);
      return;
      }
    }
  vtkErrorMacro("SetParameterScale: no parameter named " << name);
}

//----------------------------------------------------------------------------
void vtkPowellMinimizer::SetParameterScale(int i, double scale)
{
  if (i < 0 || i > this->NumberOfParameters)
    {
    vtkErrorMacro("SetParameterScale: parameter number out of range: " << i);
    return;
    }

  if (this->ParameterScales[i] != scale)
    {
    this->ParameterScales[i] = scale;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
// reset the number of parameters to zero
void vtkPowellMinimizer::Initialize()
{
  if (this->ParameterNames)
    {
    for (int i = 0; i < this->NumberOfParameters; i++)
      {
      if (this->ParameterNames[i])
        {
        delete [] this->ParameterNames[i];
        }
      }
    delete [] this->ParameterNames;
    this->ParameterNames = 0;
    }
  if (this->ParameterValues)
    {
    delete [] this->ParameterValues;
    this->ParameterValues = 0;
    }
  if (this->ParameterScales)
    {
    delete [] this->ParameterScales;
    this->ParameterScales = 0;
    }

  this->NumberOfParameters = 0;
  this->Iterations = 0;
  this->FunctionEvaluations = 0;

  this->Modified();
}

//----------------------------------------------------------------------------
void vtkPowellMinimizer::EvaluateFunction()
{
  if (this->Function)
    {
    this->Function(this->FunctionArg);
    }
  this->FunctionEvaluations++;
}

//----------------------------------------------------------------------------
double vtkPowellMinimizer::PowellGolden(
  const double *p0, double y0, const double *v, double *p, int n,
  double tol, bool *failed)
{
  // the golden ratio
  const double g = 0.6180339887498949;

  // ensure that a minimum is bracketed
  *failed = true;
  double a0 = y0;
  for (int i = 0; i < n; i++) { p[i] = p0[i] + v[i]; }
  this->EvaluateFunction();
  double a1 = this->FunctionValue;
  if (a1 < a0) { return a1; }
  for (int i = 0; i < n; i++) { p[i] = p0[i] - g*v[i]; }
  this->EvaluateFunction();
  double a2 = this->FunctionValue;
  if (a2 < a0) { return a2; }

  // do the golden section search, tol is ratio of final size to initial
  *failed = false;
  double d = 0.0;
  double a3 = a0;
  for (double r = 1.0; fabs(r) >= tol; r = g*r)
    {
    double dtry = d + r - g*r;
    for (int i = 0; i < n; i++) { p[i] = p0[i] + dtry*v[i]; }
    this->EvaluateFunction();
    a3 = this->FunctionValue;

    if (a0 < a3)
      {
      a1 = a3;
      r = -r;
      }
    else
      {
      a2 = a0;
      a0 = a3;
      d = dtry;
      }
    }

  for (int i = 0; i < n; i++) { p[i] = p0[i] + d*v[i]; }
  return a0;
}

//----------------------------------------------------------------------------
double vtkPowellMinimizer::PowellBrent(
  const double *p0, double y0, const double *vec, double *point, int n,
  double tol, bool *failed)
{
  // the golden ratio
  const double g = 0.6180339887498949;
  const double cg = 1.0 - g;

  // the bracket
  double a = -g;
  double b = 1.0;

  // ensure that a minimum is bracketed
  *failed = true;
  for (int i = 0; i < n; i++) { point[i] = p0[i] + a*vec[i]; }
  this->EvaluateFunction();
  double fa = this->FunctionValue;
  if (fa < y0) { return fa; }
  for (int i = 0; i < n; i++) { point[i] = p0[i] + b*vec[i]; }
  this->EvaluateFunction();
  double fb = this->FunctionValue;
  if (fb < y0) { return fb; }

  // else set "failed" to false (i.e. success)
  *failed = false;

  double x = 0;
  double w = 0;
  double v = 0;

  double fx = y0;
  double fw = y0;
  double fv = y0;

  double e = 0;
  double d = 1.0;

  // limit the number of line search iterations
  for (int ii = 0; ii < this->MaxIterations; ii++)
    {
    // add a fractional component to the tolerance
    double tol1 = tol + fabs(x)*1e-8;

    // midpoint
    double xc = 0.5*(a + b);
    if (fabs(x - xc) < (2*tol1 - 0.5*(b - a)))
      {
      // desired tolerance achieved
      break;
      }
    if (fabs(e) <= tol1)
      {
      // golden section
      e = ((x < xc) ? b : a) - x;
      d = cg*e;
      }
    else
      {
      // parabolic calculation
      double r = (x - w)*(fx - fv);
      double q = (x - v)*(fx - fw);
      double p = (x - v)*q + (x - w)*r;
      q = 2*(q - r);
      p = ((q > 0) ? -p : p);
      q = fabs(q);
      double etmp = e;
      e = d;

      if (p > q*(a - x) && p < q*(b - x) && fabs(p) < fabs(0.5*q*etmp))
        {
        // if we are here, the parabolic step is useful
        d = p/q;
        double u = x + d;
        if ((u - a) < 2*tol1 || (b - u) < 2*tol1)
          {
          d = ((xc - x < 0) ? -tol1 : tol1);
          }
        }
      else
        {
        // parabolic step not useful, do golden section
        e = ((x < xc) ? b : a) - x;
        d = cg*e;
        }
      }

    // compute the new position
    double u = x + d;

    // make sure step is at least as large as tol1erance
    if (fabs(d) < tol1)
      {
      // update by tol1erance
      u = x + ((d < 0) ? -tol1 : tol1);
      }

    // perform function evaluation
    for (int i = 0; i < n; i++)
      {
      point[i] = p0[i] + u*vec[i];
      }
    this->EvaluateFunction();
    double fu = this->FunctionValue;

    if (fu > fx)
      {
      if (u < x) { a = u; }
      else { b = u; }
      if (fu <= fw || w == x)
        {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
        }
      else if (fu <= fv || v == x || v == w)
        {
        v = u;
        fv = fu;
        }
      }
    else
      {
      if (u >= x) { a = x; }
      else { b = x; }
      v = w;
      w = x;
      x = u;
      fv = fw;
      fw = fx;
      fx = fu;
      }
    }

  for (int i = 0; i < n; i++)
    {
    point[i] = p0[i] + x*vec[i];
    }
  return fx;
}

//----------------------------------------------------------------------------
void vtkPowellMinimizer::PowellInitialize()
{
  int n = this->NumberOfParameters;
  double *pw = this->ParameterScales;
  delete [] this->PowellVectors;
  delete [] this->PowellWorkspace;

  // allocate memory for the current point and for
  // the conjugate directions
  double **vecs = new double *[n];
  double *work = new double[n*(n+2)];
  for (int k = 0; k < n; k++)
    {
    double *v = work + n*(k + 2);
    vecs[k] = v;
    for (int i = 0; i < n; i++) { v[i] = 0.0; }
    v[k] = pw[k];
    }

  this->PowellWorkspace = work;
  this->PowellVectors = vecs;
  this->EvaluateFunction();
}

//----------------------------------------------------------------------------
int vtkPowellMinimizer::PowellIterate()
{
  double ftol = this->Tolerance;
  double ptol = this->ParameterTolerance;
  double y = this->FunctionValue;
  double **vecs = this->PowellVectors;
  int n = this->NumberOfParameters;
  double *p = this->ParameterValues;
  double *vs = this->ParameterScales;
  double *p0 = this->PowellWorkspace;
  double *p00 = p0 + n;

  // save the current point
  for (int i = 0; i < n; i++) { p00[i] = p[i]; }
  double y00 = y;

  // go through all of the directions,
  // find the one that causes the greatest decrease
  double dymax = 0.0;
  int dymaxi = 0;
  for (int j = 0; j < n; j++)
    {
    for (int i = 0; i < n; i++) { p0[i] = p[i]; }
    double *v = vecs[j];
    bool gfailed = false;
    double y0 = y;
    // compute length of vector
    double l = 0.0;
    for (int i = 0; i < n; i++)
      {
      double w = v[i]/vs[i];
      l += w*w;
      }
    l = sqrt(l);
    double gtol = ptol/l;
    y = this->PowellBrent(p0, y0, v, p, n, gtol, &gfailed);
    double dy = y0 - y;
    if (dy > dymax)
      {
      dymax = dy;
      dymaxi = j;
      }
    }

  // compute the max distance for tolerance check
  double maxw = 0.0;
  for (int i = 0; i < n; i++)
    {
    double w = fabs(p00[i] - p[i])/vs[i];
    maxw = (maxw > w ? maxw : w);
    }

  if (2*fabs(y00 - y) <= ftol*(fabs(y00) + fabs(y)) &&
      maxw < ptol)
    {
    return 0;
    }

  // extrapolate the new point
  for (int i = 0; i < n; i++)
    {
    double ptmp = p[i];
    p[i] = 2*p[i] - p0[i];
    p0[i] = ptmp;
    }

  this->EvaluateFunction();
  double y0 = this->FunctionValue;

  // if extrapolated point is an improvement
  if (y0 < y00)
    {
    // see Numerical Recipes rationale for this in the section titled
    // "Discarding the Direction of Largest Decrease".
    double sq1 = y00 - y - dymax;
    double sq2 = y00 - y0;
    if (2*(y00 - 2*y + y0)*sq1*sq1 < sq2*sq2*dymax)
      {
      // compute the new direction
      double *vecp = vecs[dymaxi];
      for (int i = 0; i < n; i++)
        {
        vecp[i] = p0[i] - p00[i];
        }

      // move the new direction to the beginning
      vecs[dymaxi] = vecs[0];
      vecs[0] = vecp;
      }
    }

  // restore p
  for (int i = 0; i < n; i++) { p[i] = p0[i]; }
  this->FunctionValue = y;

  return 1;
}

//----------------------------------------------------------------------------
int vtkPowellMinimizer::Iterate()
{
  if (this->Iterations == 0)
    {
    if (!this->Function)
      {
      vtkErrorMacro("Iterate: Function is NULL");
      return 0;
      }
    this->PowellInitialize();
    }

  int stillgood = this->PowellIterate();
  this->Iterations++;

  return stillgood;
}

//----------------------------------------------------------------------------
void vtkPowellMinimizer::Minimize()
{
  if (this->Iterations == 0)
    {
    if (!this->Function)
      {
      vtkErrorMacro("Minimize: Function is NULL");
      return;
      }
    this->PowellInitialize();
    }

  for (; this->Iterations < this->MaxIterations; this->Iterations++)
    {
    int stillgood = this->PowellIterate();
    if (!stillgood)
      {
      break;
      }
    }
}
