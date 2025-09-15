/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImageDataStatistics.cxx,v $
  Language:  C++
  Date:      $Date: 2007/08/24 20:02:25 $
  Version:   $Revision: 1.6 $

Copyright (c) 1993-2000 Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
   of any contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 * Modified source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#include "vtkImageDataStatistics.h"

#include "vtkImageData.h"
#include "vtkObjectFactory.h"
#include "vtkCommand.h"

#if (VTK_MAJOR_VERSION >= 5)
#include "vtkInformation.h"
#include "vtkExecutive.h"
#endif

//----------------------------------------------------------------------------
vtkImageDataStatistics* vtkImageDataStatistics::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkImageDataStatistics");
  if(ret)
  {
    return (vtkImageDataStatistics*)ret;
  }
  // If the factory was unable to create the object, then create it here.
  return new vtkImageDataStatistics;
}

//----------------------------------------------------------------------------
vtkImageDataStatistics::vtkImageDataStatistics()
{
  this->AverageMagnitude = 0.0;
  this->StandardDeviation = 0.0;
  this->Count = 0;
}

//----------------------------------------------------------------------------
vtkImageDataStatistics::~vtkImageDataStatistics()
{
}

//----------------------------------------------------------------------------
void vtkImageDataStatistics::SetInput(vtkImageData *input)
{
#if (VTK_MAJOR_VERSION < 5)
  this->vtkProcessObject::SetNthInput(0, input);
#else
  // Ask the superclass to connect the input.
  this->SetNthInputConnection(0, 0, (input ? input->GetProducerPort() : 0));
#endif
}

//----------------------------------------------------------------------------
vtkImageData *vtkImageDataStatistics::GetInput()
{
#if (VTK_MAJOR_VERSION < 5)
  if (this->GetNumberOfInputs() < 1)
  {
    return NULL;
  }
  return (vtkImageData *)(this->Inputs[0]);
#else
  if (this->GetNumberOfInputConnections(0) < 1)
  {
    return NULL;
  }
  return vtkImageData::SafeDownCast(this->GetExecutive()->GetInputData(0, 0));
#endif
}

//----------------------------------------------------------------------------
template <class T>
static void vtkImageDataStatisticsExecute(vtkImageDataStatistics *self,
					  T *inPtr,
					  vtkImageData *inData,
					  double *AverageMagnitude,
					  double *StandardDeviation,
					  long int *Count)
{
  int inIdxX, inIdxY, inIdxZ;
  vtkIdType inIncX, inIncY, inIncZ;
  int wholeInExt[6];

  double sum=0.0;
  double sum_squared=0.0;
  double value;

  T *curPtr;

  // Get increments to march through data
  inData->GetWholeExtent(wholeInExt);
  inData->GetContinuousIncrements(wholeInExt, inIncX, inIncY, inIncZ);
  *Count = 0;
  // Loop through input dataset once to gather stats
  curPtr = inPtr;
  for (inIdxZ = wholeInExt[4]; inIdxZ <= wholeInExt[5]; inIdxZ++)
  {
    for (inIdxY = wholeInExt[2]; !self->AbortExecute && inIdxY <= wholeInExt[3]; inIdxY++)
    {
      for (inIdxX = wholeInExt[0]; inIdxX <= wholeInExt[1]; inIdxX++)
      {
	if (*curPtr)
 {
	  value = (double)*curPtr;
	  sum += value;
	  sum_squared += value*value;
	  (*Count)++;
 }
	else
 {
	  curPtr++;
 }
      }
      curPtr += inIncY;
    }
    curPtr += inIncZ;
  }
  *AverageMagnitude = sum / (double)*Count;
  *StandardDeviation = sqrt((sum_squared * (double)*Count - sum*sum) /
			    ((double)*Count * ((double)*Count-1.0)));
}

//----------------------------------------------------------------------------
// Description:
// Make sure input is available then call the templated execute method to
// deal with the particular data type.
void vtkImageDataStatistics::Update()
{
  vtkImageData *input = this->GetInput();
  void *inPtr;
  int wholeInExt[6];

  // make sure input is available
  if (!input)
  {
      vtkErrorMacro(<< "No input...can't execute!");
      return;
  }

  input->Update();
  input->GetWholeExtent(wholeInExt);
  inPtr = input->GetScalarPointerForExtent(wholeInExt);

  // this filter requires the input to have 1 components
  if (input->GetNumberOfScalarComponents() != 1)
  {
      vtkErrorMacro(<< "Update: input does not have 1 scalar component");
      return;
  }

  if (input->GetMTime() > this->ExecuteTime ||
      this->GetMTime() > this->ExecuteTime )
  {
      this->InvokeEvent(vtkCommand::StartEvent, NULL);

      // reset Abort flag
      this->AbortExecute = 0;
      this->Progress = 0.0;
      switch (input->GetScalarType())
      {
#if (VTK_MAJOR_VERSION < 5)
	vtkTemplateMacro6(vtkImageDataStatisticsExecute,
			  this, (VTK_TT *) (inPtr), input,
			  &(this->AverageMagnitude),
			  &(this->StandardDeviation),
			  &(this->Count));
#else
	vtkTemplateMacro(
	  vtkImageDataStatisticsExecute(this, (VTK_TT *) (inPtr), input,
					&(this->AverageMagnitude),
					&(this->StandardDeviation),
					&(this->Count)));
#endif
	default:
	  vtkErrorMacro(<< "Update: Unknown ScalarType");
	  return;
      }
      this->ExecuteTime.Modified();
      if (!this->AbortExecute)
      {
	this->UpdateProgress(1.0);
      }
      this->InvokeEvent(vtkCommand::EndEvent, NULL);
  }
  if (input->ShouldIReleaseData())
  {
    input->ReleaseData();
  }
}

//----------------------------------------------------------------------------
void vtkImageDataStatistics::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkProcessObject::PrintSelf(os,indent);

  if (!this->GetInput())
  {
    return;
  }
  os << indent << "AverageMagnitude: " << this->AverageMagnitude << "\n";
  os << indent << "StandardDeviation: " << this->StandardDeviation << "\n";
  os << indent << "Count: " << this->Count << "\n";
}

