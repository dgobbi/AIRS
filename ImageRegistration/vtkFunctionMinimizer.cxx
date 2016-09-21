/*=========================================================================

  Module: vtkFunctionMinimizer.cxx

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkFunctionMinimizer.h"

//----------------------------------------------------------------------------
vtkFunctionMinimizer::vtkFunctionMinimizer()
{
  this->Function = NULL;
  this->FunctionArg = NULL;
  this->FunctionArgDelete = NULL;

  this->NumberOfParameters = 0;
  this->ParameterNames = NULL;
  this->ParameterValues = NULL;
  this->ParameterScales = NULL;

  this->FunctionValue = 0.0;

  this->Tolerance = 1e-4;
  this->ParameterTolerance = 1e-4;
  this->MaxIterations = 1000;
  this->Iterations = 0;
  this->FunctionEvaluations = 0;
  this->AbortFlag = 0;
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
}

//----------------------------------------------------------------------------
void vtkFunctionMinimizer::PrintSelf(ostream& os, vtkIndent indent)
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
  os << indent << "AbortFlag: " << this->GetAbortFlag() << "\n";
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
double vtkFunctionMinimizer::GetParameterValue(const char *name)
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
void vtkFunctionMinimizer::SetParameterValue(const char *name, double val)
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
void vtkFunctionMinimizer::SetParameterValue(int i, double val)
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
double vtkFunctionMinimizer::GetParameterScale(const char *name)
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
void vtkFunctionMinimizer::SetParameterScale(const char *name, double scale)
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
void vtkFunctionMinimizer::SetParameterScale(int i, double scale)
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
void vtkFunctionMinimizer::Initialize()
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
  this->AbortFlag = 0;

  this->Modified();
}

//----------------------------------------------------------------------------
void vtkFunctionMinimizer::EvaluateFunction()
{
  if (!this->AbortFlag)
    {
    if (this->Function)
      {
      this->Function(this->FunctionArg);
      }
    this->FunctionEvaluations++;
    }
}

//----------------------------------------------------------------------------
int vtkFunctionMinimizer::Iterate()
{
  if (this->Iterations == 0)
    {
    if (!this->Function)
      {
      vtkErrorMacro("Iterate: Function is NULL");
      return 0;
      }
    this->Start();
    }

  if (this->AbortFlag)
    {
    return 0;
    }

  int stillgood = this->Step();

  if (this->AbortFlag)
    {
    return 0;
    }

  this->Iterations++;

  return stillgood;
}

//----------------------------------------------------------------------------
void vtkFunctionMinimizer::Minimize()
{
  if (this->Iterations == 0)
    {
    if (!this->Function)
      {
      vtkErrorMacro("Minimize: Function is NULL");
      return;
      }
    this->Start();
    }

  for (; this->Iterations < this->MaxIterations; this->Iterations++)
    {
    int stillgood = this->Step();
    if (!stillgood)
      {
      break;
      }
    }
}
