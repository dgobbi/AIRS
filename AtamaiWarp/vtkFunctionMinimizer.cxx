/*=========================================================================

  Program:   AtamaiRegistration for VTK
  Module:    $RCSfile: vtkFunctionMinimizer.cxx,v $
  Language:  C++
  Date:      $Date: 2007/08/24 20:02:25 $
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

#include "vtkFunctionMinimizer.h"
#include "vtkObjectFactory.h"

//----------------------------------------------------------------------------
static double amotry(double **p, double *y, double *ptry, double *psum, 
		     int ndim, void (*funk)(void *data), void *data,
		     double *result, int ihi, double fac)
{
  int j;
  double fac1,fac2,ytry;

  fac1 = (1.0-fac)/(double)ndim;
  fac2 = fac1-fac;
  for (j = 0; j < ndim; j++)
    {
    ptry[j] = psum[j]*fac1 - p[ihi][j]*fac2;
    }
  (*funk)(data);
  ytry = *result; 
  if (ytry < y[ihi]) 
    {
    y[ihi] = ytry;
    for (j = 0; j < ndim; j++)
      {
      psum[j] += ptry[j]-p[ihi][j];
      p[ihi][j] = ptry[j];
      }
    }
  return ytry;
}

//----------------------------------------------------------------------------
static int amoeba(double **p, double *y, double *ptry, int ndim, double ftol,
	   void (*funk)(void *data), void *data, double *result,
	   int *nfunk, int maxnfunk)
{
  int i,ihi,ilo,inhi,j,mpts;
  double rtol,sum,swap,ysave,ytry;
  double *psum = new double[ndim];
 
  mpts = ndim+1;
  *nfunk = 0;

  for (j = 0; j < ndim; j++)
    {
    sum = 0.0;
    for (i = 0; i < mpts; i++)
     {
     sum += p[i][j];
     }
    psum[j] = sum;
    }

  for (;;) 
    {
    ilo = 0;
    if (y[0] > y[1]) 
      {
      ihi = 0;
      inhi = 1;
      }
    else 
      {
      ihi = 1;
      inhi = 0;
      }
    for (i = 0; i < mpts; i++)
      {
      if (y[i] <= y[ilo]) 
	ilo = i;
      if (y[i] > y[ihi]) 
	{
	inhi = ihi;
	ihi = i;
	}
      else if (y[i] > y[inhi] && i != ihi)
	inhi = i;
      }
    
    if (fabs(y[ihi])+fabs(y[ilo]) < ftol)
      rtol = double(2.0*fabs(y[ihi]-y[ilo]));
    else
      rtol = double(2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])));

    if (rtol < ftol) {
	swap = y[1];
	y[1] = y[ilo];
	y[ilo] = swap;
	for (i = 0; i < ndim; i++) {
	  swap = p[0][i];
	  p[0][i] = p[ilo][i];
	  p[ilo][i] = swap;
	}
	break;
    }
    if (*nfunk >= maxnfunk) {
      delete [] psum;
      return -1;      /* break if greater than max number of func evals */
    }

    *nfunk += 2;
    ytry = amotry(p,y,ptry,psum,ndim,funk,data,result,ihi,double(-1.0));
    if (ytry <= y[ilo])
      ytry = amotry(p,y,ptry,psum,ndim,funk,data,result,ihi,double(2.0));
    else if (ytry >= y[inhi]) {
      ysave = y[ihi];
      ytry = amotry(p,y,ptry,psum,ndim,funk,data,result,ihi,double(0.5));
      if (ytry >= ysave) {
	for (i = 0; i < mpts; i++) {
	  if (i != ilo) {
	    for (j = 0; j < ndim; j++)
	      p[i][j] = ptry[j] = psum[j] = (p[i][j] + p[ilo][j])/double(2.0);
	    (*funk)(data);
	    y[i] = *result;
	  }
	}
	*nfunk += ndim;
	
	for (j = 0; j < ndim; j++) {
	  sum = 0;
	  for (i = 0; i < mpts; i++)
	    sum += p[i][j];
	  psum[j] = sum;
	}
      }
    }
    else
      --(*nfunk);
  }
  
  delete [] psum;
  
  return 0;
}

//----------------------------------------------------------------------------
static double minimize(double *parameters, double **vertices, int ndim, 
		       void (*funk)(void *data), void *data, double *result,
		       double tolerance, int maxiterations, int *iterations)
{
  double *y = new double[ndim+1];

  for (int k = 0; k < ndim+1; k++) 
    {
    for (int l = 0; l < ndim; l++)
      {
      parameters[l] = vertices[k][l];
      }

    (*funk)(data);
    y[k] = *result;
    }

  amoeba(vertices,y,parameters,ndim,tolerance,funk,data,result,
	 iterations,maxiterations);
  *result = y[1]; // copy the lowest result in the *result
  delete [] y;

  //set x equal to lowest of amoeba vertices
  for (int j = 0; j < ndim; j++)
    parameters[j] = vertices[0][j];

  return *result;
}

//----------------------------------------------------------------------------
void vtkFunctionMinimizerFunction(void *data)
{
  vtkFunctionMinimizer *self = (vtkFunctionMinimizer *)data;
  /* old code for when a vtkFunctionParser was used
  for (int i = 0; i < self->NumberOfParameters; i++)
    {
    self->Function->SetScalarVariableValue(self->ParameterIndices[i],
    }
  self->ScalarResult = self->Function->GetScalarResult();
  */
  if (self->Function)
    {
    self->Function(self->FunctionArg);
    }
}					   

//----------------------------------------------------------------------------
vtkFunctionMinimizer *vtkFunctionMinimizer::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkFunctionMinimizer");
  if(ret)
    {
    return (vtkFunctionMinimizer*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkFunctionMinimizer;
}

//----------------------------------------------------------------------------
vtkFunctionMinimizer::vtkFunctionMinimizer()
{
  this->Function = NULL;
  this->FunctionArg = NULL;
  this->FunctionArgDelete = NULL;

  this->NumberOfParameters = 0;
  this->ParameterNames = NULL;
  this->ParameterIndices = NULL;
  this->Parameters = NULL;
  this->ParameterBrackets = NULL;
  this->Vertices = NULL;

  this->ScalarResult = 0.0;

  this->Tolerance = 0.005;
  this->MaxIterations = 1000;
  this->Iterations = 0;
}
  
//----------------------------------------------------------------------------
vtkFunctionMinimizer::~vtkFunctionMinimizer()
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
  if (this->ParameterIndices)
    {
    delete [] this->ParameterIndices;
    this->ParameterIndices = NULL;
    }
  if (this->Parameters)
    {
    delete [] this->Parameters;
    this->Parameters = NULL;
    }
  if (this->ParameterBrackets)
    {
    delete [] this->ParameterBrackets;
    this->ParameterBrackets = NULL;
    }
  if (this->Vertices)
    {
    delete [] *this->Vertices;
    delete [] this->Vertices;
    this->Vertices = NULL;
    }

  this->NumberOfParameters = 0;
} 

//----------------------------------------------------------------------------
void vtkFunctionMinimizer::PrintSelf(ostream& os, vtkIndent indent)
{
  this->vtkObject::PrintSelf(os, indent);
  os << indent << "ScalarResult: " << this->ScalarResult << "\n";
  os << indent << "MaxIterations: " << this->MaxIterations << "\n";
  os << indent << "Iterations: " << this->Iterations << "\n";
  os << indent << "Tolerance: " << this->Tolerance << "\n";
}

//----------------------------------------------------------------------------
void vtkFunctionMinimizer::SetFunction(void (*f)(void *), void *arg)
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
void vtkFunctionMinimizer::SetFunctionArgDelete(void (*f)(void *))
{
  if ( f != this->FunctionArgDelete)
    {
    this->FunctionArgDelete = f;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
double *vtkFunctionMinimizer::GetScalarVariableBracket(const char *name)
{
  static double errval[2] = { 0.0, 0.0 };

  for (int i = 0; i < this->NumberOfParameters; i++)
    {
    if (strcmp(name,this->ParameterNames[i]) == 0)
      {
      return this->ParameterBrackets[i];
      }
    }

  vtkErrorMacro("GetScalarVariableBracket: no parameter named " << name);
  return errval;
}

//----------------------------------------------------------------------------
double *vtkFunctionMinimizer::GetScalarVarPtr()
{
  return this->Parameters;
}

//----------------------------------------------------------------------------
double vtkFunctionMinimizer::GetScalarVariableValue(const char *name)
{
  for (int i = 0; i < this->NumberOfParameters; i++)
    {
    if (strcmp(name,this->ParameterNames[i]) == 0)
      {
      return this->Parameters[i];
      }
    }
  vtkErrorMacro("GetScalarVariableValue: no parameter named " << name);
  return 0.0;
}

//----------------------------------------------------------------------------
// initialize the simplex, also find the indices of the variables
int vtkFunctionMinimizer::Initialize()
{
  if (!this->Function)
    {
    vtkErrorMacro("Initialize: Funtion is NULL");
    return 0;
    }

  for (int l = 0; l < this->NumberOfParameters; l++)
    {
      // initial parameter values are middle of bracket
      this->Vertices[0][l] = 0.5*(this->ParameterBrackets[l][0] +
				  this->ParameterBrackets[l][1]);

      // set up the simplex vertices
      for (int m = 1; m <= this->NumberOfParameters; m++)
	{
	  this->Vertices[m][l] = this->Vertices[0][l];
	  if ((m-1)==l) this->Vertices[m][l]=this->ParameterBrackets[l][1];
	}
    }
  
  this->Iterations = 0;
  
  return 1;
}

//----------------------------------------------------------------------------
void vtkFunctionMinimizer::SetScalarVariableBracket(const char *name, 
						    double bmin, double bmax)
{
  int i;

  for (i = 0; i < this->NumberOfParameters; i++)
    {
    if (strcmp(name,this->ParameterNames[i]) == 0)
      {
      if (this->ParameterBrackets[i][0] != bmin ||
	  this->ParameterBrackets[i][1] != bmax)
	{
	this->ParameterBrackets[i][0] = bmin;
	this->ParameterBrackets[i][1] = bmax;
//  	if ((!finite(bmin)) | (!finite(bmax))) {
//  	  cout << "dying in SetScalarVariableBracket\n";
//  	  int *seg =0;
//  	  *seg = 3;}
	  
	this->Modified();
	}
      return;
      }
    }

  int n = this->NumberOfParameters + 1;
  char **newParameterNames = new char *[n];
  int *newParameterIndices = new int[n];
  double *newParameters = new double[n];
  double (*newParameterBrackets)[2] = new double[n][2];
  double **newVertices = new double *[n+1];
  double *mem = new double[n*(n+1)];

  for (i = 0; i < this->NumberOfParameters; i++)
    {
    newParameterNames[i] = this->ParameterNames[i];
    newParameterBrackets[i][0] = this->ParameterBrackets[i][0];
    newParameterBrackets[i][1] = this->ParameterBrackets[i][1];
    newVertices[i] = &mem[n*i];
    }

  char *cp = new char[strlen(name)+8];
  strcpy(cp,name);
  newParameterNames[n-1] = cp;
  newParameterBrackets[n-1][0] = bmin;
  newParameterBrackets[n-1][1] = bmax;
  newVertices[n-1] = &mem[n*(n-1)];
  newVertices[n] = &mem[n*n];

  if (this->ParameterNames)
    {
    delete [] this->ParameterNames;
    }
  if (this->ParameterIndices)
    {
    delete [] this->ParameterIndices;
    }
  if (this->Parameters)
    {
    delete [] this->Parameters;
    }
  if (this->ParameterBrackets)
    {
    delete [] this->ParameterBrackets;
    }
  if (this->Vertices)
    {
    delete [] *this->Vertices;
    delete [] this->Vertices;
    }

  this->NumberOfParameters = n;
  this->ParameterNames = newParameterNames;
  this->ParameterIndices = newParameterIndices;
  this->Parameters = newParameters;
  this->ParameterBrackets = newParameterBrackets;
  this->Vertices = newVertices;

  this->Modified();
}

//----------------------------------------------------------------------------
void vtkFunctionMinimizer::Minimize()
{
  if (!this->Initialize())
    {
    return;
    }
  minimize(this->Parameters, this->Vertices, this->NumberOfParameters,
	   &vtkFunctionMinimizerFunction, this, &this->ScalarResult,
	   this->Tolerance, this->MaxIterations, &this->Iterations);
}

