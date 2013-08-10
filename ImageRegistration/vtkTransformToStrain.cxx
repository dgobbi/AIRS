/*=========================================================================

  Program:   Atamai Image Registration and Segmentation
  Module:    vtkTransformToStrain.cxx

  Copyright (c) 2013 David Gobbi
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

  * Neither the name of David Gobbi, nor the names of any authors nor
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=========================================================================*/
#include "vtkTransformToStrain.h"

#include "vtkMath.h"
#include "vtkAbstractTransform.h"
#include "vtkIdentityTransform.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"

vtkStandardNewMacro(vtkTransformToStrain);

vtkCxxSetObjectMacro(vtkTransformToStrain,Input,vtkAbstractTransform);

//----------------------------------------------------------------------------
vtkTransformToStrain::vtkTransformToStrain()
{
  this->Input = NULL;

  this->OutputScalarType = VTK_FLOAT;
  this->OutputValue = GreensStrain;

  for (int i = 0; i < 3; i++)
    {
    this->OutputExtent[2*i] = this->OutputExtent[2*i+1] = 0;
    this->OutputOrigin[i] = 0.0;
    this->OutputSpacing[i] = 1.0;
    }

  this->ValueScale = 1.0;
  this->ValueShift = 0.0;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkTransformToStrain::~vtkTransformToStrain()
{
  this->SetInput(static_cast<vtkAbstractTransform*>(0));
}

//----------------------------------------------------------------------------
void vtkTransformToStrain::PrintSelf(ostream& os, vtkIndent indent)
{
  int i;

  this->Superclass::PrintSelf(os,indent);

  os << indent << "Input: (" << this->Input << ")\n";

  os << indent << "OutputSpacing: (" << this->OutputSpacing[0];
  for (i = 1; i < 3; ++i)
    {
    os << ", " << this->OutputSpacing[i];
    }
  os << ")\n";

  os << indent << "OutputOrigin: (" << this->OutputOrigin[0];
  for (i = 1; i < 3; ++i)
    {
    os << ", " << this->OutputOrigin[i];
    }
  os << ")\n";

  os << indent << "OutputExtent: (" << this->OutputExtent[0];
  for (i = 1; i < 6; ++i)
    {
    os << ", " << this->OutputExtent[i];
    }
  os << ")\n";

  os << indent << "OutputScalarType: " <<
    vtkImageScalarTypeNameMacro(this->OutputScalarType) << "\n";

  os << indent << "OutputValue: " << this->GetOutputValueAsString() << "\n";

  os << indent << "ValueScale: " << this->ValueScale << "\n";
  os << indent << "ValueShift: " << this->ValueShift << "\n";
}

//----------------------------------------------------------------------------
const char *vtkTransformToStrain::GetOutputValueAsString()
{
  const char *v = "";

  switch (this->OutputValue)
    {
    case GreensStrain:
      v = "GreensStrain";
      break;
    case DeformationGradient:
      v = "DeformationGradient";
      break;
    }

  return v;
}

//----------------------------------------------------------------------------
void vtkTransformToStrain::ComputeGreensStrain(
  const double F[3][3], double G[3][3])
{
  vtkMath::Transpose3x3(F, G);
  vtkMath::Multiply3x3(G, F, G);

  G[0][0] -= 1.0;
  G[1][1] -= 1.0;
  G[2][2] -= 1.0;

  G[0][0] *= 0.5;
  G[0][1] *= 0.5;
  G[0][2] *= 0.5;
  G[1][0] *= 0.5;
  G[1][1] *= 0.5;
  G[1][2] *= 0.5;
  G[2][0] *= 0.5;
  G[2][1] *= 0.5;
  G[2][2] *= 0.5;
}

//----------------------------------------------------------------------------
// This method returns the largest data that can be generated.
void vtkTransformToStrain::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  if (this->GetInput() == NULL)
    {
    vtkErrorMacro("Missing input");
    return;
    }

  // update the transform, maybe in the future make transforms part of the
  // pipeline
  this->Input->Update();

  outInfo->Set(
    vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), this->OutputExtent, 6);

  outInfo->Set(vtkDataObject::SPACING(), this->OutputSpacing, 3);
  outInfo->Set(vtkDataObject::ORIGIN(), this->OutputOrigin, 3);

  vtkDataObject::SetPointDataActiveScalarInfo(
    outInfo, this->OutputScalarType, 9);
}

//----------------------------------------------------------------------------
// Return the maximum and minimum value of the tensor elements over
// the entire output extent -- this is extremely robust and extremely
// inefficient, it should be possible to do much better than this.
void vtkTransformToStrainMinMax(
  vtkTransformToStrain *self, int extent[6], int operation,
  double &minValue, double &maxValue)
{
  vtkAbstractTransform *transform = self->GetInput();
  transform->Update();

  if (!transform)
    {
    minValue = -1.0;
    maxValue = +1.0;
    return;
    }

  double *spacing = self->GetOutputSpacing();
  double *origin = self->GetOutputOrigin();

  maxValue = -1e37;
  minValue = +1e37;

  double point[3], newPoint[3], tensor[3][3];

  for (int k = extent[4]; k <= extent[5]; k++)
    {
    point[2] = k*spacing[2] + origin[2];
    for (int j = extent[2]; j <= extent[3]; j++)
      {
      point[1] = j*spacing[1] + origin[1];
      for (int i = extent[0]; i <= extent[1]; i++)
        {
        point[0] = i*spacing[0] + origin[0];

        transform->InternalTransformDerivative(
          point, newPoint, tensor);

        if (operation == vtkTransformToStrain::GreensStrain)
          {
          vtkTransformToStrain::ComputeGreensStrain(tensor, tensor);
          }

        for (int l = 0; l < 3; l++)
          {
          for (int m = 0; m < 3; m++)
            {
            double v = tensor[l][m];

            if (v > maxValue)
              {
              maxValue = v;
              }

            if (v < minValue)
              {
              minValue = v;
              }
            }
          }
        }
      }
    }
}

//----------------------------------------------------------------------------
void vtkTransformToStrain::UpdateShiftScale()
{
  int outputType = this->OutputScalarType;

  // nothing to do for double or double
  if (outputType == VTK_DOUBLE || outputType == VTK_FLOAT)
    {
    this->ValueShift = 0.0;
    this->ValueScale = 1.0;
    vtkDebugMacro(<< "value (scale, shift) = (" <<
                  this->ValueScale << ", " <<
                  this->ValueShift << ")");
    return;
    }

  // check mtime
  if (this->ShiftScaleTime.GetMTime() > this->GetMTime())
    {
    return;
    }

  // get the maximum displacement
  int operation = this->OutputValue;
  double minDisplacement, maxDisplacement;
  vtkTransformToStrainMinMax(
    this, this->OutputExtent, operation, minDisplacement, maxDisplacement);

  vtkDebugMacro(<< "displacement (min, max) = (" <<
                minDisplacement << ", " << maxDisplacement << ")");

  double typeMin,typeMax;

  switch (outputType)
    {
    case VTK_SHORT:
      typeMin = VTK_SHORT_MIN;
      typeMax = VTK_SHORT_MAX;
      break;
    case VTK_UNSIGNED_SHORT:
      typeMin = VTK_UNSIGNED_SHORT_MIN;
      typeMax = VTK_UNSIGNED_SHORT_MAX;
      break;
    case VTK_CHAR:
      typeMin = VTK_CHAR_MIN;
      typeMax = VTK_CHAR_MAX;
      break;
    case VTK_UNSIGNED_CHAR:
      typeMin = VTK_UNSIGNED_CHAR_MIN;
      typeMax = VTK_UNSIGNED_CHAR_MAX;
      break;
    default:
      vtkErrorMacro(<< "UpdateShiftScale: Unknown input ScalarType");
      return;
    }

  this->ValueScale = ((maxDisplacement - minDisplacement)/
                             (typeMax - typeMin));
  this->ValueShift = ((typeMax*minDisplacement-typeMin*maxDisplacement)/
                             (typeMax - typeMin));

  if (this->ValueScale == 0.0)
    {
    this->ValueScale = 1.0;
    }

  vtkDebugMacro(<< "displacement (scale, shift) = (" <<
                this->ValueScale << ", " <<
                this->ValueShift << ")");

  this->ShiftScaleTime.Modified();
}

//----------------------------------------------------------------------------
// macros to ensure proper round-to-nearest behaviour

inline void vtkOutputRound(double val, unsigned char& rnd)
{
  val = (val < -128.0 ? -128.0 : val);
  val = (val > 127.0 ? 127.0 : val);
  int ival = static_cast<int>(val + 128.5) - 128;
  rnd = static_cast<signed char>(ival);
}

inline void vtkOutputRound(double val, signed char& rnd)
{
  val = (val < 0 ? 0 : val);
  val = (val > 255.0 ? 255.0 : val);
  rnd = static_cast<unsigned char>(val + 0.5);
}

inline void vtkOutputRound(double val, short& rnd)
{
  val = (val < -32768.0 ? -32768.0 : val);
  val = (val > 32767.0 ? 32767.0 : val);
  int ival = static_cast<int>(val + 32768.5) - 32768;
  rnd = static_cast<short>(ival);
}

inline void vtkOutputRound(double val, unsigned short& rnd)
{
  val = (val < 0 ? 0 : val);
  val = (val > 65535.0 ? 65535.0 : val);
  rnd = static_cast<unsigned short>(val + 0.5);
}

inline void vtkOutputRound(double val, double& rnd)
{
  rnd = static_cast<double>(val);
}

inline void vtkOutputRound(double val, float& rnd)
{
  rnd = static_cast<float>(val);
}

//----------------------------------------------------------------------------
template<class T>
void vtkTransformToStrainExecute(
  vtkTransformToStrain *self, vtkImageData *output, T *outPtr, int extent[6],
  double shift, double scale, int operation, int id)
{
  vtkAbstractTransform *transform = self->GetInput();
  int isIdentity = 0;
  if (transform == 0)
    {
    transform = vtkIdentityTransform::New();
    isIdentity = 1;
    }

  double *spacing = output->GetSpacing();
  double *origin = output->GetOrigin();
  vtkIdType *increments = output->GetIncrements();

  double invScale = 1.0/scale;

  double point[3];
  double newPoint[3];
  double tensor[3][3];

  T *outPtr0 = outPtr;

  unsigned long count = 0;
  unsigned long target = static_cast<unsigned long>(
    (extent[5] - extent[4] + 1)*(extent[3] - extent[2] + 1)/50.0);
  target++;

  for (int k = extent[4]; k <= extent[5]; k++)
    {
    point[2] = k*spacing[2] + origin[2];
    T *outPtr1 = outPtr0;

    for (int j = extent[2]; j <= extent[3]; j++)
      {
      if (id == 0)
        {
        if (count % target == 0)
          {
          self->UpdateProgress(count/(50.0*target));
          }
        count++;
        }

      point[1] = j*spacing[1] + origin[1];
      outPtr = outPtr1;

      for (int i = extent[0]; i <= extent[1]; i++)
        {
        point[0] = i*spacing[0] + origin[0];

        transform->InternalTransformDerivative(
          point, newPoint, tensor);

        if (operation == vtkTransformToStrain::GreensStrain)
          {
          vtkTransformToStrain::ComputeGreensStrain(tensor, tensor);
          }

        for (int ii = 0; ii < 3; ii++)
          {
          for (int jj = 0; jj < 3; jj++)
            {
            vtkOutputRound((tensor[ii][jj] - shift)*invScale, *outPtr++);
            }
          }
        }

      outPtr1 += increments[1];
      }

    outPtr0 += increments[2];
    }

  if (isIdentity)
    {
    transform->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkTransformToStrain::RequestData(
  vtkInformation* vtkNotUsed( request ),
  vtkInformationVector** vtkNotUsed( inputVector ),
  vtkInformationVector* outputVector)
{
  // get the data object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkImageData *output = vtkImageData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  output->SetExtent(
    outInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()));
  output->AllocateScalars();
  int *extent = output->GetExtent();

  void *outPtr = output->GetScalarPointerForExtent(extent);
  int outputType = output->GetScalarType();
  int operation = this->OutputValue;

  this->UpdateShiftScale();

  double scale = this->ValueScale;
  double shift = this->ValueShift;

  int id = 0;

  switch (outputType)
    {
    case VTK_FLOAT:
      vtkTransformToStrainExecute(
        this, output, static_cast<float *>(outPtr), extent,
        shift, scale, operation, id);
      break;
    case VTK_DOUBLE:
      vtkTransformToStrainExecute(
        this, output, static_cast<double *>(outPtr), extent,
        shift, scale, operation, id);
      break;
    case VTK_SHORT:
      vtkTransformToStrainExecute(
        this, output, static_cast<short *>(outPtr), extent,
        shift, scale, operation, id);
      break;
    case VTK_UNSIGNED_SHORT:
      vtkTransformToStrainExecute(
        this, output, static_cast<unsigned short *>(outPtr), extent,
        shift, scale, operation, id);
      break;
    case VTK_SIGNED_CHAR:
      vtkTransformToStrainExecute(
        this, output, static_cast<unsigned char *>(outPtr), extent,
        shift, scale, operation, id);
      break;
    case VTK_UNSIGNED_CHAR:
      vtkTransformToStrainExecute(
        this, output, static_cast<unsigned char *>(outPtr), extent,
        shift, scale, operation, id);
      break;
    default:
      vtkErrorMacro(<< "Execute: Unknown output ScalarType");
    }
}

//----------------------------------------------------------------------------
unsigned long vtkTransformToStrain::GetMTime()
{
  unsigned long mtime = this->Superclass::GetMTime();

  if (this->Input)
    {
    unsigned long mtime2 = this->Input->GetMTime();
    if (mtime2 > mtime)
      {
      mtime = mtime2;
      }
    }

  return mtime;
}

//----------------------------------------------------------------------------
int vtkTransformToStrain::ProcessRequest(
  vtkInformation* request, vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // generate the data
  if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
    {
    this->RequestData(request, inputVector, outputVector);
    return 1;
    }

  // execute information
  if(request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
    {
    this->RequestInformation(request, inputVector, outputVector);
    // after executing set the origin and spacing from the info
    for (int i = 0; i < this->GetNumberOfOutputPorts(); ++i)
      {
      vtkInformation* info = outputVector->GetInformationObject(i);
      vtkImageData *output =
        vtkImageData::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));
      // if execute info didn't set origin and spacing then we set them
      if (!info->Has(vtkDataObject::ORIGIN()))
        {
        info->Set(vtkDataObject::ORIGIN(), 0, 0, 0);
        info->Set(vtkDataObject::SPACING(), 1, 1, 1);
        }
      if (output)
        {
        output->SetOrigin(info->Get(vtkDataObject::ORIGIN()));
        output->SetSpacing(info->Get(vtkDataObject::SPACING()));
        }
      }
    return 1;
    }

  return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

//----------------------------------------------------------------------------
vtkImageData* vtkTransformToStrain::GetOutput()
{
  return vtkImageData::SafeDownCast(this->GetOutputDataObject(0));
}

int vtkTransformToStrain::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // now add our info
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
  return 1;
}
