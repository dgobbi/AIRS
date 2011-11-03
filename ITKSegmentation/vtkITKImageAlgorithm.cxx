/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkITKImageAlgorithm.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkITKImageAlgorithm.h"

#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkObjectFactory.h"



vtkCxxRevisionMacro(vtkITKImageAlgorithm, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkITKImageAlgorithm);

//----------------------------------------------------------------------------
vtkITKImageAlgorithm::vtkITKImageAlgorithm()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfInputPorts(1);
}
vtkITKImageAlgorithm::~vtkITKImageAlgorithm()
{
}

//----------------------------------------------------------------------------
void vtkITKImageAlgorithm::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkITKImageAlgorithm::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector))
{
  int inExt[6];
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkImageData *inputData =
    vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  inputData->GetWholeExtent(inExt);
  inputData->SetUpdateExtent(inExt);

  return 1;
}

//----------------------------------------------------------------------------
int vtkITKImageAlgorithm::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *vtkNotUsed(outputVector))
{
  return 1;
}
