/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkNelderMeadMinimizer.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkNelderMeadMinimizer.h"
#include "vtkObjectFactory.h"

#define  N_STEPS_NO_VALUE_IMPROVEMENT  2
#define  N_STEPS_NO_PARAM_IMPROVEMENT  18

vtkStandardNewMacro(vtkNelderMeadMinimizer);

//----------------------------------------------------------------------------
vtkNelderMeadMinimizer::vtkNelderMeadMinimizer()
{
  this->ContractionRatio = 0.5;
  this->ExpansionRatio = 2.0;
  this->AmoebaVertices = NULL;
  this->AmoebaValues = NULL;
  this->AmoebaSum = NULL;
  this->AmoebaSize = 0;
  this->AmoebaHighValue = 0;
  this->AmoebaNStepsNoImprovement = 0;
}

//----------------------------------------------------------------------------
vtkNelderMeadMinimizer::~vtkNelderMeadMinimizer()
{
  this->TerminateAmoeba();
}

//----------------------------------------------------------------------------
void vtkNelderMeadMinimizer::PrintSelf(ostream& os, vtkIndent indent)
{
  os << indent << "ContractionRatio: " << this->GetContractionRatio() << "\n";
  os << indent << "ExpansionRatio: " << this->GetExpansionRatio() << "\n";
}

//----------------------------------------------------------------------------
int vtkNelderMeadMinimizer::CheckParameterTolerance()
{
  int n = this->NumberOfParameters;

  double *vertex0 = this->AmoebaVertices[0];
  double *scales = this->ParameterScales;
  double size = 0;

  for (int i = 1; i <= n; i++)
  {
    double *vertex = this->AmoebaVertices[i];
    for (int j = 0; j < n; j++)
    {
      double d = fabs((vertex[j] - vertex0[j])/scales[j]);
      size = ((d < size) ? size : d);
    }
  }

  if (size != this->AmoebaSize)
  {
    this->AmoebaNStepsNoImprovement = N_STEPS_NO_VALUE_IMPROVEMENT-1;
  }
  this->AmoebaSize = size;
  // if amoeba is static, only make a set number of tries
  if (this->AmoebaNStepsNoImprovement >
      (N_STEPS_NO_VALUE_IMPROVEMENT + N_STEPS_NO_PARAM_IMPROVEMENT))
  {
    return 1;
  }

  return (size <= this->ParameterTolerance);
}

//----------------------------------------------------------------------------
void vtkNelderMeadMinimizer::Start()
{
  this->AmoebaSize = 0;
  this->InitializeAmoeba();
}

//----------------------------------------------------------------------------
int vtkNelderMeadMinimizer::Step()
{
  int improved = this->PerformAmoeba();
  int paramsWithinTol = 0;
  if (!improved)
  {
    paramsWithinTol = this->CheckParameterTolerance();
  }
  this->GetAmoebaParameterValues();

  return (improved || !paramsWithinTol);
}

/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

/* ------------------------------------------------
  This code has been modified from the original.  Several macros
  have been expanded, functions have been renamed to match VTK
  conventions, and the formatting has been changed.
*/

/* ----------------------------- MNI Header -----------------------------------
@NAME       : vtkAmoebaNumericallyClose
@INPUT      : n1
              n2
              threshold_ratio
@OUTPUT     :
@RETURNS    : true if the numbers are within the threshold ratio
@DESCRIPTION: Decides if two numbers are close to each other.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    :         1993    David MacDonald
@MODIFIED   :         2002    David Gobbi
---------------------------------------------------------------------------- */

#define  VTK_AMOEBA_SMALLEST  1.0e-20

static  int  vtkAmoebaNumericallyClose(double  n1,
                                       double  n2,
                                       double  threshold_ratio )
{
  double  avg, diff, abs_n1, abs_n2;

  diff = n1 - n2;
  if( diff < 0.0 )
  {
    diff = -diff;
  }

  abs_n1 = (n1 < 0.0 ? -n1 : n1);
  abs_n2 = (n2 < 0.0 ? -n2 : n2);

  if( abs_n1 < VTK_AMOEBA_SMALLEST || abs_n2 < VTK_AMOEBA_SMALLEST )
  {
    return( abs_n1 < threshold_ratio && abs_n2 < threshold_ratio );
  }

  avg = (n1 + n2) / 2.0;

  if( avg == 0.0 )
  {
    return( diff <= threshold_ratio );
  }

  if( avg < 0.0 )
  {
    avg = -avg;
  }

  return( (diff / avg) <= threshold_ratio );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : InitializeAmoeba
@INPUT      :
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Initializes the amoeba structure to minimize the function.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    :         1993    David MacDonald
@MODIFIED   :         2002    David Gobbi
---------------------------------------------------------------------------- */

void  vtkNelderMeadMinimizer::InitializeAmoeba()
{
  int    i, j;

  this->TerminateAmoeba();

  int n_parameters = this->NumberOfParameters;
  this->AmoebaNStepsNoImprovement = 0;
  this->AmoebaVertices = new double *[n_parameters+1];
  this->AmoebaVertices[0] = new double[n_parameters*(n_parameters+1)];

  for( i = 1 ; i < n_parameters+1 ; i++)
  {
    this->AmoebaVertices[i] = this->AmoebaVertices[i-1] + n_parameters;
  }

  this->AmoebaValues = new double[n_parameters+1];

  this->AmoebaSum = new double[n_parameters];

  for (j = 0; j < n_parameters; j++)
  {
    this->AmoebaSum[j] = 0.0;
  }

  for( i = 0 ; i < n_parameters+1 ; i++ )
  {
    for( j = 0; j < n_parameters ; j++ )
    {
      this->AmoebaVertices[i][j] = this->ParameterValues[j];
      if( i > 0 && j == i - 1 )
      {
        this->AmoebaVertices[i][j] =
          this->ParameterValues[j] + this->ParameterScales[j];
      }
      this->AmoebaSum[j] += this->ParameterValues[j];
    }
  }
  for( i = 0 ; i < n_parameters+1 ; i++ )
  {
    for( j = 0; j < n_parameters; j++ )
    {
      this->ParameterValues[j] = this->AmoebaVertices[i][j];
    }
    this->EvaluateFunction();
    this->AmoebaValues[i] = this->FunctionValue;
  }

  for ( j = 0 ; j < n_parameters ; j++ )
  {
    this->ParameterValues[j] = this->AmoebaVertices[0][j];
  }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GetAmoebaParameterValues
@INPUT      :
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Passes back the current position of the amoeba (best value),
              and returns the function value at that point.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    :         1993    David MacDonald
@MODIFIED   :         2002    David Gobbi
---------------------------------------------------------------------------- */

void vtkNelderMeadMinimizer::GetAmoebaParameterValues()
{
  int   i, j, low;

  low = 0;
  for( i = 1 ; i < this->NumberOfParameters+1 ; i++ )
  {
    if( this->AmoebaValues[i] < this->AmoebaValues[low] )
    {
      low = i;
    }
  }

  for( j = 0 ; j < this->NumberOfParameters ; j++ )
  {
    this->ParameterValues[j] = this->AmoebaVertices[low][j];
  }

  this->FunctionValue = this->AmoebaValues[low];
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : TerminateAmoeba
@INPUT      :
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Frees the amoeba.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    :         1993    David MacDonald
@MODIFIED   :         2002    David Gobbi
---------------------------------------------------------------------------- */

void  vtkNelderMeadMinimizer::TerminateAmoeba()
{
  if (this->AmoebaVertices)
  {
    delete [] this->AmoebaVertices[0];
    delete [] this->AmoebaVertices;
    this->AmoebaVertices = NULL;
  }
  if (this->AmoebaValues)
  {
    delete [] this->AmoebaValues;
    this->AmoebaValues = NULL;
  }
  if (this->AmoebaSum)
  {
    delete [] this->AmoebaSum;
    this->AmoebaSum = NULL;
  }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : TryAmoeba
@INPUT      : sum
              high
              fac
@OUTPUT     :
@RETURNS    : value
@DESCRIPTION: Does a modification to the high vertex of the amoeba and
              returns the value of the new point.  If the new point is
              better (smaller value), it replaces the high vertex of the
              amoeba.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    :         1993    David MacDonald
@MODIFIED   :         2002    David Gobbi
---------------------------------------------------------------------------- */

double  vtkNelderMeadMinimizer::TryAmoeba(double  sum[],
                                      int     high,
                                      double  fac )
{
  int    j;
  double y_try, fac1, fac2;
  double  *parameters;

  parameters = this->ParameterValues;

  fac1 = (1.0 - fac) / this->NumberOfParameters;
  fac2 = fac - fac1;

  for( j = 0 ; j < this->NumberOfParameters ; j++ )
  {
    parameters[j] = (sum[j] * fac1 + this->AmoebaVertices[high][j] * fac2);
  }

  this->EvaluateFunction();
  y_try = this->FunctionValue;

  if( y_try < this->AmoebaValues[high] )
  {
    this->AmoebaValues[high] = y_try;
    for( j = 0 ; j < this->NumberOfParameters ; j++ )
    {
      sum[j] += parameters[j] - this->AmoebaVertices[high][j];
      this->AmoebaVertices[high][j] = parameters[j];
    }
  }

  return( y_try );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : PerformAmoeba
@INPUT      :
@OUTPUT     :

@RETURNS    : true if numerically significant improvement
@DESCRIPTION: Performs one iteration of an amoeba, returning true if a
              numerically significant improvement has been found recently.
              Even if it returns 0, you can keep calling this function,
              since it may be contracting with no improvement, but will
              eventually shrink small enough to get an improvment.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    :         1993    David MacDonald
@MODIFIED   :         2002    David Gobbi
---------------------------------------------------------------------------- */

int vtkNelderMeadMinimizer::PerformAmoeba()
{
  int      i, j, low, high, next_high;
  double   y_try, y_save;
  int      improvement_found;

  improvement_found = 1;

  if( this->AmoebaValues[0] > this->AmoebaValues[1] )
  {
    high = 0;
    next_high = 1;
  }
  else
  {
    high = 1;
    next_high = 0;
  }

  low = next_high;

  for( i = 2 ; i < this->NumberOfParameters+1 ; i++ )
  {
    if( this->AmoebaValues[i] < this->AmoebaValues[low] )
    {
      low = i;
    }
    else if( this->AmoebaValues[i] > this->AmoebaValues[high] )
    {
      next_high = high;
      high = i;
    }
    else if( this->AmoebaValues[i] > this->AmoebaValues[next_high] )
    {
      next_high = i;
    }
  }

  if( this->AmoebaValues[high] == this->AmoebaHighValue ||
      vtkAmoebaNumericallyClose( this->AmoebaValues[low],
                                 this->AmoebaValues[high],
                                 this->Tolerance ) )
  {
    ++this->AmoebaNStepsNoImprovement;
    if( this->AmoebaNStepsNoImprovement >= N_STEPS_NO_VALUE_IMPROVEMENT )
    {
      improvement_found = 0;
    }
  }
  else
  {
    this->AmoebaNStepsNoImprovement = 0;
  }

  this->AmoebaHighValue = this->AmoebaValues[high];

  y_try = this->TryAmoeba( this->AmoebaSum, high, -1.0 );

  if( y_try <= this->AmoebaValues[low] )
  {
    TryAmoeba( this->AmoebaSum, high, this->ExpansionRatio );
  }
  else if( y_try >= this->AmoebaValues[next_high] )
  {
    y_save = this->AmoebaValues[high];
    y_try = TryAmoeba( this->AmoebaSum, high, this->ContractionRatio );

    if( y_try >= y_save )
    {
      for( i = 0 ; i < this->NumberOfParameters+1 ; i++)
      {
        if( i != low )
        {
          for( j = 0 ; j < this->NumberOfParameters ; j++ )
          {
            this->ParameterValues[j] = (this->AmoebaVertices[i][j] +
                                        this->AmoebaVertices[low][j]) / 2.0;
            this->AmoebaVertices[i][j] = this->ParameterValues[j];
          }

          this->EvaluateFunction();
          this->AmoebaValues[i] = this->FunctionValue;
        }
      }

      for( j = 0 ; j < this->NumberOfParameters ; j++ )
      {
        this->AmoebaSum[j] = 0.0;
        for( i = 0 ; i < this->NumberOfParameters+1 ; i++ )
        {
          this->AmoebaSum[j] += this->AmoebaVertices[i][j];
        }
      }
    }
  }

  return( improvement_found );
}
