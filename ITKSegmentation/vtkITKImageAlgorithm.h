/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkITKImageAlgorithm.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkITKImageAlgorithm - superclass for VTK image algorithms that use ITK internally.
// .SECTION Description

#ifndef __vtkITKImageAlgorithm_h
#define __vtkITKImageAlgorithm_h


#include "vtkImageAlgorithm.h"


class vtkITKImageAlgorithm : public vtkImageAlgorithm
{
public:
  static vtkITKImageAlgorithm *New();
  vtkTypeRevisionMacro(vtkITKImageAlgorithm, vtkImageAlgorithm);

  virtual void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkITKImageAlgorithm();
  ~vtkITKImageAlgorithm();

  virtual int RequestData(vtkInformation *,
			  vtkInformationVector **,
			  vtkInformationVector *);
  virtual int RequestUpdateExtent(vtkInformation *, vtkInformationVector **,
                                  vtkInformationVector *);

private:
  vtkITKImageAlgorithm(const vtkITKImageAlgorithm&);  // Not implemented.
  void operator=(const vtkITKImageAlgorithm&);  // Not implemented.
};

#endif
