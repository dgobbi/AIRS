/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSimpleShapeStencilSource.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSimpleShapeStencilSource.h"

#include "vtkMath.h"
#include "vtkDataArray.h"
#include "vtkImageData.h"
#include "vtkImageStencilData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <math.h>

vtkStandardNewMacro(vtkSimpleShapeStencilSource);
vtkCxxSetObjectMacro(vtkSimpleShapeStencilSource, InformationInput, vtkImageData);

//----------------------------------------------------------------------------
vtkSimpleShapeStencilSource::vtkSimpleShapeStencilSource()
{
  this->SetNumberOfInputPorts(0);

  this->Shape = vtkSimpleShapeStencilSource::BOX;

  this->Center[0] = 0.0;
  this->Center[1] = 0.0;
  this->Center[2] = 0.0;

  this->Size[0] = 1.0;
  this->Size[1] = 1.0;
  this->Size[2] = 1.0;

  this->InformationInput = NULL;

  this->OutputOrigin[0] = 0;
  this->OutputOrigin[1] = 0;
  this->OutputOrigin[2] = 0;

  this->OutputSpacing[0] = 1;
  this->OutputSpacing[1] = 1;
  this->OutputSpacing[2] = 1;

  this->OutputWholeExtent[0] = 0;
  this->OutputWholeExtent[1] = 0;
  this->OutputWholeExtent[2] = 0;
  this->OutputWholeExtent[3] = 0;
  this->OutputWholeExtent[4] = 0;
  this->OutputWholeExtent[5] = 0;
}

//----------------------------------------------------------------------------
vtkSimpleShapeStencilSource::~vtkSimpleShapeStencilSource()
{
  this->SetInformationInput(NULL);
}

//----------------------------------------------------------------------------
void vtkSimpleShapeStencilSource::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "InformationInput: " << this->InformationInput << "\n";

  os << indent << "OutputSpacing: " << this->OutputSpacing[0] << " " <<
    this->OutputSpacing[1] << " " << this->OutputSpacing[2] << "\n";
  os << indent << "OutputOrigin: " << this->OutputOrigin[0] << " " <<
    this->OutputOrigin[1] << " " << this->OutputOrigin[2] << "\n";
  os << indent << "OutputWholeExtent: " << this->OutputWholeExtent[0] << " " <<
    this->OutputWholeExtent[1] << " " << this->OutputWholeExtent[2] << " " <<
    this->OutputWholeExtent[3] << " " << this->OutputWholeExtent[4] << " " <<
    this->OutputWholeExtent[5] << "\n";

  os << indent << "Shape: " << this->GetShapeAsString() << "\n";
  os << indent << "Center: " << this->Center[0] << " "
     << this->Center[1] << " " << this->Center[2] << "\n";
  os << indent << "Size: " << this->Size[0] << " "
     << this->Size[1] << " " << this->Size[2] << "\n";
}

//----------------------------------------------------------------------------
const char *vtkSimpleShapeStencilSource::GetShapeAsString()
{
  switch (this->Shape)
    {
    case vtkSimpleShapeStencilSource::BOX:
      return "Box";
    case vtkSimpleShapeStencilSource::ELLIPSOID:
      return "Ellipsoid";
    case vtkSimpleShapeStencilSource::CYLINDERX:
      return "CylinderX";
    case vtkSimpleShapeStencilSource::CYLINDERY:
      return "CylinderY";
    case vtkSimpleShapeStencilSource::CYLINDERZ:
      return "CylinderZ";
    }
  return "";
}

//----------------------------------------------------------------------------
// tolerance for stencil operations

#define VTK_STENCIL_TOL 7.62939453125e-06

//----------------------------------------------------------------------------
static int vtkSimpleShapeStencilSourceBox(
  vtkSimpleShapeStencilSource *self, vtkImageStencilData *data,
  int extent[6], double origin[3], double spacing[3])
{
  // for keeping track of progress
  unsigned long count = 0;
  unsigned long target = static_cast<unsigned long>(
    (extent[5] - extent[4] + 1)*(extent[3] - extent[2] + 1)/50.0);
  target++;

  double center[3];
  self->GetCenter(center);

  double size[3];
  self->GetSize(size);

  for (int i = 0; i < 3; i++)
    {
    double r = 0.5*size[i];
    if (r < 0) { r = -r; }
    double dmin = (center[i] - r - origin[i])/spacing[i];
    double dmax = (center[i] + r - origin[i])/spacing[i];
    int emin = vtkMath::Floor(dmin + VTK_STENCIL_TOL);
    int emax = -vtkMath::Floor(-dmax + VTK_STENCIL_TOL);
    if (extent[2*i] < emin)
      {
      extent[2*i] = emin;
      }
    if (extent[2*i+1] > emax)
      {
      extent[2*i+1] = emax;
      }
    }

  for (int idZ = extent[4]; idZ <= extent[5]; idZ++)
    {
    for (int idY = extent[2]; idY <= extent[3]; idY++)
      {
      if (count%target == 0) 
        {
        self->UpdateProgress(count/(50.0*target));
        }
      count++;

      int r1 = extent[0];
      int r2 = extent[1];

      if (r2 >= r1)
        {
        data->InsertNextExtent(r1, r2, idY, idZ);
        }
      } // for idY
    } // for idZ

  return 1;
}

//----------------------------------------------------------------------------
static int vtkSimpleShapeStencilSourceEllipsoid(
  vtkSimpleShapeStencilSource *self, vtkImageStencilData *data,
  int extent[6], double origin[3], double spacing[3])
{
  // for keeping track of progress
  unsigned long count = 0;
  unsigned long target = static_cast<unsigned long>(
    (extent[5] - extent[4] + 1)*(extent[3] - extent[2] + 1)/50.0);
  target++;

  double center[3];
  self->GetCenter(center);

  double size[3];
  self->GetSize(size);

  double icenter[3];
  double iradius[3];

  for (int i = 0; i < 3; i++)
    {
    icenter[i] = (center[i] - origin[i])/spacing[i];
    iradius[i] = 0.5*size[i]/spacing[i];
    if (iradius[i] < 0) { iradius[i] = -iradius[i]; }
    int emin = vtkMath::Floor(icenter[i] - iradius[i] + VTK_STENCIL_TOL);
    int emax = -vtkMath::Floor(-icenter[i] - iradius[i] + VTK_STENCIL_TOL);
    if (extent[2*i] < emin)
      {
      extent[2*i] = emin;
      }
    if (extent[2*i+1] > emax)
      {
      extent[2*i+1] = emax;
      }
    }

  for (int idZ = extent[4]; idZ <= extent[5]; idZ++)
    {
    double z = (idZ - icenter[2])/iradius[2];

    for (int idY = extent[2]; idY <= extent[3]; idY++)
      {
      if (count%target == 0) 
        {
        self->UpdateProgress(count/(50.0*target));
        }
      count++;

      double y = (idY - icenter[1])/iradius[1];
      double d = 1.0 - y*y - z*z;
      if (d < -VTK_STENCIL_TOL)
        {
        continue;
        }
      else if (d < 0)
        {
        d = 0;
        }

      double x = sqrt(d);
      int r1 = vtkMath::Floor(icenter[0] - x*iradius[0] + VTK_STENCIL_TOL);
      int r2 = -vtkMath::Floor(-icenter[0] - x*iradius[0] + VTK_STENCIL_TOL);

      if (r1 < extent[0])
        {
        r1 = extent[0];
        }
      if (r2 > extent[1])
        {
        r2 = extent[1];
        }

      if (r2 >= r1)
        {
        data->InsertNextExtent(r1, r2, idY, idZ);
        }
      } // for idY
    } // for idZ

  return 1;
}

//----------------------------------------------------------------------------
static int vtkSimpleShapeStencilSourceCylinderX(
  vtkSimpleShapeStencilSource *self, vtkImageStencilData *data,
  int extent[6], double origin[3], double spacing[3])
{
  // for keeping track of progress
  unsigned long count = 0;
  unsigned long target = static_cast<unsigned long>(
    (extent[5] - extent[4] + 1)*(extent[3] - extent[2] + 1)/50.0);
  target++;

  double center[3];
  self->GetCenter(center);

  double size[3];
  self->GetSize(size);

  double icenter[3];
  double iradius[3];

  for (int i = 0; i < 3; i++)
    {
    icenter[i] = (center[i] - origin[i])/spacing[i];
    iradius[i] = 0.5*size[i]/spacing[i];
    if (iradius[i] < 0) { iradius[i] = -iradius[i]; }
    int emin = vtkMath::Floor(icenter[i] - iradius[i] + VTK_STENCIL_TOL);
    int emax = -vtkMath::Floor(-icenter[i] - iradius[i] + VTK_STENCIL_TOL);
    if (extent[2*i] < emin)
      {
      extent[2*i] = emin;
      }
    if (extent[2*i+1] > emax)
      {
      extent[2*i+1] = emax;
      }
    }

  for (int idZ = extent[4]; idZ <= extent[5]; idZ++)
    {
    double z = (idZ - icenter[2])/iradius[2];

    for (int idY = extent[2]; idY <= extent[3]; idY++)
      {
      if (count%target == 0) 
        {
        self->UpdateProgress(count/(50.0*target));
        }
      count++;

      double y = (idY - icenter[1])/iradius[1];
      if (y*y + z*z > 1.0 + VTK_STENCIL_TOL)
        {
        continue;
        }

      int r1 = extent[0];
      int r2 = extent[1];

      if (r2 >= r1)
        {
        data->InsertNextExtent(r1, r2, idY, idZ);
        }
      } // for idY
    } // for idZ

  return 1;
}

//----------------------------------------------------------------------------
static int vtkSimpleShapeStencilSourceCylinderY(
  vtkSimpleShapeStencilSource *self, vtkImageStencilData *data,
  int extent[6], double origin[3], double spacing[3])
{
  // for keeping track of progress
  unsigned long count = 0;
  unsigned long target = static_cast<unsigned long>(
    (extent[5] - extent[4] + 1)*(extent[3] - extent[2] + 1)/50.0);
  target++;

  double center[3];
  self->GetCenter(center);

  double size[3];
  self->GetSize(size);

  double icenter[3];
  double iradius[3];

  for (int i = 0; i < 3; i++)
    {
    icenter[i] = (center[i] - origin[i])/spacing[i];
    iradius[i] = 0.5*size[i]/spacing[i];
    if (iradius[i] < 0) { iradius[i] = -iradius[i]; }
    int emin = vtkMath::Floor(icenter[i] - iradius[i] + VTK_STENCIL_TOL);
    int emax = -vtkMath::Floor(-icenter[i] - iradius[i] + VTK_STENCIL_TOL);
    if (extent[2*i] < emin)
      {
      extent[2*i] = emin;
      }
    if (extent[2*i+1] > emax)
      {
      extent[2*i+1] = emax;
      }
    }

  for (int idZ = extent[4]; idZ <= extent[5]; idZ++)
    {
    double z = (idZ - icenter[2])/iradius[2];

    for (int idY = extent[2]; idY <= extent[3]; idY++)
      {
      if (count%target == 0) 
        {
        self->UpdateProgress(count/(50.0*target));
        }
      count++;

      double d = 1.0 - z*z;
      if (d < -VTK_STENCIL_TOL)
        {
        continue;
        }
      else if (d < 0)
        {
        d = 0;
        }

      double x = sqrt(d);
      int r1 = vtkMath::Floor(icenter[0] - x*iradius[0] + VTK_STENCIL_TOL);
      int r2 = -vtkMath::Floor(-icenter[0] - x*iradius[0] + VTK_STENCIL_TOL);

      if (r1 < extent[0])
        {
        r1 = extent[0];
        }
      if (r2 > extent[1])
        {
        r2 = extent[1];
        }

      if (r2 >= r1)
        {
        data->InsertNextExtent(r1, r2, idY, idZ);
        }
      } // for idY
    } // for idZ

  return 1;
}

//----------------------------------------------------------------------------
static int vtkSimpleShapeStencilSourceCylinderZ(
  vtkSimpleShapeStencilSource *self, vtkImageStencilData *data,
  int extent[6], double origin[3], double spacing[3])
{
  // for keeping track of progress
  unsigned long count = 0;
  unsigned long target = static_cast<unsigned long>(
    (extent[5] - extent[4] + 1)*(extent[3] - extent[2] + 1)/50.0);
  target++;

  double center[3];
  self->GetCenter(center);

  double size[3];
  self->GetSize(size);

  double icenter[3];
  double iradius[3];

  for (int i = 0; i < 3; i++)
    {
    icenter[i] = (center[i] - origin[i])/spacing[i];
    iradius[i] = 0.5*size[i]/spacing[i];
    if (iradius[i] < 0) { iradius[i] = -iradius[i]; }
    int emin = vtkMath::Floor(icenter[i] - iradius[i] + VTK_STENCIL_TOL);
    int emax = -vtkMath::Floor(-icenter[i] - iradius[i] + VTK_STENCIL_TOL);
    if (extent[2*i] < emin)
      {
      extent[2*i] = emin;
      }
    if (extent[2*i+1] > emax)
      {
      extent[2*i+1] = emax;
      }
    }

  for (int idZ = extent[4]; idZ <= extent[5]; idZ++)
    {
    for (int idY = extent[2]; idY <= extent[3]; idY++)
      {
      if (count%target == 0) 
        {
        self->UpdateProgress(count/(50.0*target));
        }
      count++;

      double y = (idY - icenter[1])/iradius[1];
      double d = 1.0 - y*y;
      if (d < -VTK_STENCIL_TOL)
        {
        continue;
        }
      else if (d < 0)
        {
        d = 0;
        }

      double x = sqrt(d);
      int r1 = vtkMath::Floor(icenter[0] - x*iradius[0] + VTK_STENCIL_TOL);
      int r2 = -vtkMath::Floor(-icenter[0] - x*iradius[0] + VTK_STENCIL_TOL);

      if (r1 < extent[0])
        {
        r1 = extent[0];
        }
      if (r2 > extent[1])
        {
        r2 = extent[1];
        }

      if (r2 >= r1)
        {
        data->InsertNextExtent(r1, r2, idY, idZ);
        }
      } // for idY
    } // for idZ

  return 1;
}


//----------------------------------------------------------------------------
int vtkSimpleShapeStencilSource::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  int extent[6];
  double origin[3];
  double spacing[3];
  int result = 1;

  this->Superclass::RequestData(request, inputVector, outputVector);

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkImageStencilData *data = vtkImageStencilData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent);
  outInfo->Get(vtkDataObject::ORIGIN(), origin);
  outInfo->Get(vtkDataObject::SPACING(), spacing);

  switch (this->Shape)
    {
    case vtkSimpleShapeStencilSource::BOX:
      result = vtkSimpleShapeStencilSourceBox(
        this, data, extent, origin, spacing);
      break;
    case vtkSimpleShapeStencilSource::ELLIPSOID:
      result = vtkSimpleShapeStencilSourceEllipsoid(
        this, data, extent, origin, spacing);
      break;
    case vtkSimpleShapeStencilSource::CYLINDERX:
      result = vtkSimpleShapeStencilSourceCylinderX(
        this, data, extent, origin, spacing);
      break;
    case vtkSimpleShapeStencilSource::CYLINDERY:
      result = vtkSimpleShapeStencilSourceCylinderY(
        this, data, extent, origin, spacing);
      break;
    case vtkSimpleShapeStencilSource::CYLINDERZ:
      result = vtkSimpleShapeStencilSourceCylinderZ(
        this, data, extent, origin, spacing);
      break;
    }

  return result;
}

//----------------------------------------------------------------------------
int vtkSimpleShapeStencilSource::RequestInformation(
  vtkInformation *,
  vtkInformationVector **,
  vtkInformationVector *outputVector)
{
  int wholeExtent[6];
  double spacing[3];
  double origin[3];

  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  for (int i = 0; i < 3; i++)
    {
    wholeExtent[2*i] = this->OutputWholeExtent[2*i];
    wholeExtent[2*i+1] = this->OutputWholeExtent[2*i+1];
    spacing[i] = this->OutputSpacing[i];
    origin[i] = this->OutputOrigin[i];
    }

  // If InformationInput is set, then get the spacing,
  // origin, and whole extent from it.
  if (this->InformationInput)
    {
    this->InformationInput->UpdateInformation();
    this->InformationInput->GetWholeExtent(wholeExtent);
    this->InformationInput->GetSpacing(spacing);
    this->InformationInput->GetOrigin(origin);
    }

  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
               wholeExtent, 6);
  outInfo->Set(vtkDataObject::SPACING(), spacing, 3);
  outInfo->Set(vtkDataObject::ORIGIN(), origin, 3);

  return 1;
}
