/*=========================================================================

  Module: vtkPowellMinimizer.cxx

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkPowellMinimizer.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkPowellMinimizer);

//----------------------------------------------------------------------------
vtkPowellMinimizer::vtkPowellMinimizer()
{
  this->PowellWorkspace = 0;
  this->PowellVectors = 0;
}

//----------------------------------------------------------------------------
vtkPowellMinimizer::~vtkPowellMinimizer()
{
  delete [] this->PowellVectors;
  delete [] this->PowellWorkspace;
}

//----------------------------------------------------------------------------
void vtkPowellMinimizer::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
double vtkPowellMinimizer::PowellBrent(
  const double *p0, double y0, const double *vec, double *point, int n,
  const double bracket[3], double tol)
{
  // the golden ratio
  const double g = 0.6180339887498949;
  const double cg = 1.0 - g;

  // the bracket
  double a = bracket[0];
  double b = bracket[2];
  if (a > b)
    {
    b = bracket[0];
    a = bracket[2];
    }

  double x = bracket[1];
  double w = x;
  double v = x;

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

    // make sure step is at least as large as tolerance
    if (fabs(d) < tol1)
      {
      // update by tolerance
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
double vtkPowellMinimizer::PowellBracket(
  const double *p0, double y0, const double *vec, double *point, int n,
  double bracket[3], bool *failed)
{
  // the golden ratio
  const double g = 0.6180339887498949;
  // maximum growth allowed
  const double growlim = 110;

  double fa = y0;
  double xa = 0.0;
  double xb = 1.0;

  for (int i = 0; i < n; i++)
    {
    point[i] = p0[i] + xb*vec[i];
    }
  this->EvaluateFunction();
  double fb = this->FunctionValue;
  if (fa < fb)
    {
    xa = xb;
    xb = 0.0;
    fa = fb;
    fb = y0;
    }

  double xc = xb + g*(xb - xa);

  for (int i = 0; i < n; i++)
    {
    point[i] = p0[i] + xc*vec[i];
    }
  this->EvaluateFunction();
  double fc = this->FunctionValue;

  int ii = 0;
  while (fc < fb)
    {
    double tmp1 = (xb - xa)*(fb - fc);
    double tmp2 = (xb - xc)*(fb - fa);
    double val = tmp2 - tmp1;
    const double tinyval = 1e-21;
    val = ((val < tinyval) ? tinyval : val);
    double denom = 2*val;
    double w = xb - ((xb - xc)*tmp2 - (xb - xa)*tmp1)/denom;
    double wlim = xb + growlim*(xc - xb);
    if (ii++ > this->MaxIterations)
      {
      break;
      }

    double fw = 0;
    if ((w - xc)*(xb - w) > 0)
      {
      for (int i = 0; i < n; i++)
        {
        point[i] = p0[i] + w*vec[i];
        }
      this->EvaluateFunction();
      fw = this->FunctionValue;
      if (fw < fc)
        {
        xa = xb;
        xb = w;
        fa = fb;
        fb = fw;
        break;
        }
      else if (fw > fb)
        {
        xc = w;
        fc = fw;
        break;
        }
      w = xc + g*(xc - xb);
      for (int i = 0; i < n; i++)
        {
        point[i] = p0[i] + w*vec[i];
        }
      this->EvaluateFunction();
      fw = this->FunctionValue;
      }
    else if ((w - wlim)*(wlim - xc) >= 0)
      {
      w = wlim;
      for (int i = 0; i < n; i++)
        {
        point[i] = p0[i] + w*vec[i];
        }
      this->EvaluateFunction();
      fw = this->FunctionValue;
      }
    else if ((w - wlim)*(xc - w) >= 0)
      {
      for (int i = 0; i < n; i++)
        {
        point[i] = p0[i] + w*vec[i];
        }
      this->EvaluateFunction();
      fw = this->FunctionValue;
      if (fw < fc)
        {
        xb = xc;
        xc = w;
        w = xc + g*(xc - xb);
        fb = fc;
        fc = fw;
        for (int i = 0; i < n; i++)
          {
          point[i] = p0[i] + w*vec[i];
          }
        this->EvaluateFunction();
        fw = this->FunctionValue;
        }
      }
    else
      {
      w = xc + g*(xc - xb);
      for (int i = 0; i < n; i++)
        {
        point[i] = p0[i] + w*vec[i];
        }
      this->EvaluateFunction();
      fw = this->FunctionValue;
      }
    xa = xb;
    xb = xc;
    xc = w;
    fa = fb;
    fb = fc;
    fc = fw;
    }

  bracket[0] = xa;
  bracket[1] = xb;
  bracket[2] = xc;

  if (fa <= fb)
    {
    *failed = true;
    xb = xa;
    fb = fa;
    }
  else if (fc <= fb)
    {
    *failed = true;
    xb = xc;
    fb = fc;
    }
  else
    {
    *failed = false;
    }

  for (int i = 0; i < n; i++)
    {
    point[i] = p0[i] + xb*vec[i];
    }
  return fb;
}

//----------------------------------------------------------------------------
void vtkPowellMinimizer::Start()
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
int vtkPowellMinimizer::Step()
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
    double bracket[3];
    bool failed = false;
    y = this->PowellBracket(p0, y0, v, p, n, bracket, &failed);
    if (!failed)
      {
      y = this->PowellBrent(p0, y, v, p, n, bracket, gtol);
      }
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
