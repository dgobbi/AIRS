/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCalcCrossCorrelation.cxx,v $
  Language:  C++
  Date:      $Date: 2007/08/24 20:02:25 $
  Version:   $Revision: 1.4 $


Copyright (c) 1993-2001 Ken Martin, Will Schroeder, Bill Lorensen 
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
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#include "vtkCalcCrossCorrelation.h"

#include "vtkObjectFactory.h"
#include "vtkCommand.h"
#include "vtkImageData.h"
#include "vtkImageStencilData.h"

#if (VTK_MAJOR_VERSION >= 5) 
#include "vtkInformation.h"
#include "vtkExecutive.h"
#endif

//---------------------------------------------------------------------------
vtkCalcCrossCorrelation* vtkCalcCrossCorrelation::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkCalcCrossCorrelation");
  if(ret)
    {
    return (vtkCalcCrossCorrelation*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkCalcCrossCorrelation;
}

// Constructs with initial 0 values.
vtkCalcCrossCorrelation::vtkCalcCrossCorrelation()
{
  this->CrossCorrelation = 0.0;
  this->ReverseStencil = 0;

#if (VTK_MAJOR_VERSION >= 5) 
  // we have the image inputs and the optional stencil input
  this->SetNumberOfInputPorts(2);
#endif
}

// Destroy any allocated memory.
vtkCalcCrossCorrelation::~vtkCalcCrossCorrelation()
{
}

// Description:
// Specifies the input datasets
void vtkCalcCrossCorrelation::SetInput1(vtkImageData *input)
{
#if (VTK_MAJOR_VERSION >= 5)
  this->SetNthInputConnection(0, 0, (input ? input->GetProducerPort() : 0));
#else
  this->vtkProcessObject::SetNthInput(0, input);
#endif
}
void vtkCalcCrossCorrelation::SetInput2(vtkImageData *input)
{
#if (VTK_MAJOR_VERSION >= 5)
  this->SetNthInputConnection(0, 1, (input ? input->GetProducerPort() : 0));
#else
  this->vtkProcessObject::SetNthInput(1, input);
#endif
}

//----------------------------------------------------------------------------
void vtkCalcCrossCorrelation::SetStencil(vtkImageStencilData *stencil)
{
#if (VTK_MAJOR_VERSION >= 5)
  // if stencil is null, then set the input port to null
  this->SetNthInputConnection(1, 0, 
    (stencil ? stencil->GetProducerPort() : 0));
#else
  this->vtkProcessObject::SetNthInput(2, stencil);
#endif
}

//----------------------------------------------------------------------------
vtkImageData *vtkCalcCrossCorrelation::GetInput1()
{
#if (VTK_MAJOR_VERSION >= 5)
  if (this->GetNumberOfInputConnections(0) < 1)
    {
    return NULL;
    }
  return vtkImageData::SafeDownCast(this->GetExecutive()->GetInputData(0, 0));
#else
  if (this->GetNumberOfInputs() < 1)    
    {
    return NULL;
    }
  return (vtkImageData *)(this->Inputs[0]);
#endif
}
//----------------------------------------------------------------------------
vtkImageData *vtkCalcCrossCorrelation::GetInput2()
{
#if (VTK_MAJOR_VERSION >= 5)
  if (this->GetNumberOfInputConnections(0) < 2)
    {
    return NULL;
    }
  return vtkImageData::SafeDownCast(this->GetExecutive()->GetInputData(0, 1));
#else
  if (this->GetNumberOfInputs() < 2)    
    {
    return NULL;
    }
  return (vtkImageData *)(this->Inputs[1]);
#endif
}

//----------------------------------------------------------------------------
vtkImageStencilData *vtkCalcCrossCorrelation::GetStencil()
{
#if (VTK_MAJOR_VERSION >= 5)
  if (this->GetNumberOfInputConnections(1) < 1)
    {
    return NULL;
    }
  return vtkImageStencilData::SafeDownCast(
    this->GetExecutive()->GetInputData(1, 0));
#else
  if (this->NumberOfInputs < 3)
    {
    return NULL;
    }
  else
    {
    return (vtkImageStencilData *)(this->Inputs[2]);
    }
#endif
}

// Description:
// Make sure input is available then call up execute method...
void vtkCalcCrossCorrelation::Update()
{
  vtkImageData *input1 = this->GetInput1();
  vtkImageData *input2 = this->GetInput2();
  
  // make sure input is available
  if ( ! input1 || ! input2)
    {
    vtkErrorMacro(<< "No inputs...can't execute!");
    return;
    }

  input1->UpdateData();
  input2->UpdateData();

  vtkImageStencilData *stencil = this->GetStencil();
  if (stencil)
    {
    stencil->SetSpacing(input1->GetSpacing());
    stencil->SetOrigin(input1->GetOrigin());
    stencil->SetUpdateExtent(stencil->GetWholeExtent());
    stencil->Update();
    }

  if (input1->GetMTime() > this->ExecuteTime || input2->GetMTime() > this->ExecuteTime || 
      this->GetMTime() > this->ExecuteTime )
    {
    if ( input1->GetDataReleased() )
      {
      input1->Update();
      }
    if ( input2->GetDataReleased() )
      {
      input2->Update();
      }
    this->InvokeEvent(vtkCommand::StartEvent,NULL);

    // reset Abort flag
    this->AbortExecute = 0;
    this->Progress = 0.0;
    this->Execute();
    this->ExecuteTime.Modified();
    if ( !this->AbortExecute )
      {
      this->UpdateProgress(1.0);
      }

    this->InvokeEvent(vtkCommand::EndEvent,NULL);
    }
  if ( input1->ShouldIReleaseData() )
    {
    input1->ReleaseData();
    }
  if ( input2->ShouldIReleaseData() )
    {
    input2->ReleaseData();
    }
}

template <class T>
void vtkCorrelationHelper(vtkCalcCrossCorrelation *self,
			  vtkImageData *inData1, T *in1Ptr, 
			  vtkImageData *inData2, T *in2Ptr,
			  double *Correlation)
{
  int idX, idY, idZ;
  vtkIdType incX, incY, incZ;
  int minX, maxX, minY, maxY, minZ, maxZ;
  int pminX, pmaxX;
  int iter;
  *Correlation = 0.0;
  double topSum = 0;
  double Sum1 = 0;
  double Sum2 = 0;
  T *temp1Ptr;
  T *temp2Ptr;

  vtkImageStencilData *stencil = self->GetStencil();

  inData1->GetUpdateExtent(minX, maxX, minY, maxY, minZ, maxZ);
  inData1->GetIncrements(incX, incY, incZ);
  
  unsigned long count = 0;
  unsigned long target = (unsigned long)((maxZ - minZ + 1)*(maxY - minY +1)/50.0);
  target++;

  for (idZ = minZ; idZ <= maxZ; idZ++)
    {
     for (idY = minY; idY <= maxY; idY++)
      {
      if (!(count%target)) 
	{
        self->UpdateProgress(count/(50.0*target));
	}
      count++;

      // loop over stencil sub-extents
      iter = 0;
      if (self->GetReverseStencil())
	{// flag that we want the complementary extents
	iter = -1;
	}

      pminX = minX;
      pmaxX = maxX;
      while ((stencil !=0 &&
	      stencil->GetNextExtent(pminX, pmaxX, minX, maxX, idY, idZ, iter)) ||
	     (stencil == 0 && iter++ == 0))
	{
	// set up pointers to the sub-extents
	temp1Ptr = in1Ptr + (incZ*(idZ - minZ) +
			     incY*(idY - minY) +
			     (pminX-minX));
	temp2Ptr = in2Ptr + (incZ*(idZ - minZ) +
			     incY*(idY - minY) +
			     (pminX-minX));
	// compute over the sub-extent
	for (idX = pminX; idX <= pmaxX; idX++)
	  { // only add pixels if both are non-zero
	  if (*temp1Ptr && *temp2Ptr) 
	    {
	    topSum += (double)*temp1Ptr * (double)*temp2Ptr;
	    Sum1   += (double)*temp1Ptr * (double)*temp1Ptr;
	    Sum2   += (double)*temp2Ptr * (double)*temp2Ptr;
	    }
	  temp1Ptr++;
	  temp2Ptr++;
	  }
	}
      }
    }
  *Correlation = (double) topSum / (sqrt((double)Sum1) * sqrt((double(Sum2)))); 
}

// Description:
// This method computes normalized cross correlation.
// Requires that the inputs are the same extent.
void vtkCalcCrossCorrelation::Execute()
{
  vtkImageData *input1 = this->GetInput1();
  vtkImageData *input2 = this->GetInput2();

  if (!input1 || !input2)
    {
    vtkErrorMacro("Must set both inputs");
    return;
    }
  if (input1->GetScalarType() != input2->GetScalarType())
    {
    vtkErrorMacro("Scalar Types must be the same");
    return;
    }

  vtkImageStencilData *stencil = this->GetStencil();

  if (stencil)
    {
    stencil->SetWholeExtent(input1->GetWholeExtent());
    stencil->SetUpdateExtent(stencil->GetWholeExtent());
    stencil->Update();
    }
   
  input1->SetUpdateExtent(input1->GetWholeExtent());
  input1->Update();
  input2->SetUpdateExtent(input2->GetWholeExtent());
  input2->Update();
  void *in1Ptr = input1->GetScalarPointer();
  void *in2Ptr = input2->GetScalarPointer();

  if (!in1Ptr || !in2Ptr)
    {
    vtkErrorMacro("No data in the ImageDatas");
    return;
    }

  switch (input1->GetScalarType())
    {
#if (VTK_MAJOR_VERSION < 5)
    vtkTemplateMacro6(vtkCorrelationHelper,this,
		      input1, (VTK_TT *)(in1Ptr),
		      input2, (VTK_TT *)(in2Ptr),
		      &this->CrossCorrelation);
#else
    vtkTemplateMacro(
      vtkCorrelationHelper(this,
			   input1, (VTK_TT *)(in1Ptr),
			   input2, (VTK_TT *)(in2Ptr),
			   &this->CrossCorrelation));
#endif
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType (must be char, us char, int, us int");
      return;
    }
}

void vtkCalcCrossCorrelation::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkProcessObject::PrintSelf(os,indent);

  os << indent << "Input 1: "<< this->GetInput1() << "\n";
  os << indent << "Input 2: "<< this->GetInput2() << "\n";
  os << indent << "Stencil: " << this->GetStencil() << "\n";
  os << indent << "ReverseStencil: " << (this->ReverseStencil ? "On\n" : "Off\n");
  os << indent << "CrossCorrelation: " << this->GetCrossCorrelation () << "\n";
}

#if (VTK_MAJOR_VERSION >= 5) 
//----------------------------------------------------------------------------
int vtkCalcCrossCorrelation::FillInputPortInformation(int port, 
                                                      vtkInformation* info)
{
  if (port == 0)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    }
  if (port == 1)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageStencilData");
    // the stencil input is optional
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }
  return 1;
}
#endif











