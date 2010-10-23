/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSimpleShapeStencilSource.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSimpleShapeStencilSource - create simple mask shapes
// .SECTION Description
// vtkSimpleShapeStencilSource will create an image stencil with a
// simple shape like a box, a sphere, or a cylinder.  Its output can
// be used with vtkImageStecil or other vtk classes that apply
// a stencil to an image.
// .SECTION see also
// vtkImplicitFunctionToImageStencil vtkPolyDataToImageStencil

#ifndef __vtkSimpleShapeStencilSource_h
#define __vtkSimpleShapeStencilSource_h


#include "vtkImageStencilSource.h"

class vtkImageData;

class VTK_EXPORT vtkSimpleShapeStencilSource : public vtkImageStencilSource
{
public:
  static vtkSimpleShapeStencilSource *New();
  vtkTypeMacro(vtkSimpleShapeStencilSource, vtkImageStencilSource);
  void PrintSelf(ostream& os, vtkIndent indent);

//BTX
  enum {
    BOX = 0,
    ELLIPSOID = 1,
    CYLINDERX = 2,
    CYLINDERY = 3,
    CYLINDERZ = 4
  };
//ETX

  // Description:
  // The shape of the region of interest.  Cylinders can be oriented
  // along the X, Y, or Z axes.
  vtkGetMacro(Shape, int);
  vtkSetClampMacro(Shape, int, BOX, CYLINDERZ);
  void SetShapeToBox() { this->SetShape(BOX); };
  void SetShapeToEllipsoid() { this->SetShape(ELLIPSOID); };
  void SetShapeToCylinderX() { this->SetShape(CYLINDERX); };
  void SetShapeToCylinderY() { this->SetShape(CYLINDERY); };
  void SetShapeToCylinderZ() { this->SetShape(CYLINDERZ); };
  virtual const char *GetShapeAsString();

  // Description:
  // The center of the region of interest.
  vtkGetVector3Macro(Center, double);
  vtkSetVector3Macro(Center, double);

  // Description:
  // The size of the region of interest.
  vtkGetVector3Macro(Size, double);
  vtkSetVector3Macro(Size, double);

  // Description:
  // Set a vtkImageData that has the Spacing, Origin, and
  // WholeExtent that will be used for the stencil.  This
  // input should be set to the image that you wish to
  // apply the stencil to.  If you use this method, then
  // any values set with the SetOutputSpacing, SetOutputOrigin,
  // and SetOutputWholeExtent methods will be ignored.
  virtual void SetInformationInput(vtkImageData*);
  vtkGetObjectMacro(InformationInput, vtkImageData);

  // Description:
  // Set the Origin to be used for the stencil.  It should be
  // set to the Origin of the image you intend to apply the
  // stencil to. The default value is (0,0,0).
  vtkSetVector3Macro(OutputOrigin, double);
  vtkGetVector3Macro(OutputOrigin, double);

  // Description:
  // Set the Spacing to be used for the stencil. It should be
  // set to the Spacing of the image you intend to apply the
  // stencil to. The default value is (1,1,1)
  vtkSetVector3Macro(OutputSpacing, double);
  vtkGetVector3Macro(OutputSpacing, double);

  // Description:
  // Set the whole extent for the stencil (anything outside
  // this extent will be considered to be "outside" the stencil).
  // If this is not set, then the stencil will always use
  // the requested UpdateExtent as the stencil extent.
  vtkSetVector6Macro(OutputWholeExtent, int);
  vtkGetVector6Macro(OutputWholeExtent, int);  

protected:
  vtkSimpleShapeStencilSource();
  ~vtkSimpleShapeStencilSource();

  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  
  int Shape;
  double Center[3];
  double Size[3];

  vtkImageData *InformationInput;

  int OutputWholeExtent[6];
  double OutputOrigin[3];
  double OutputSpacing[3];

private:
  vtkSimpleShapeStencilSource(const vtkSimpleShapeStencilSource&);  // Not implemented.
  void operator=(const vtkSimpleShapeStencilSource&);  // Not implemented.
};

#endif
