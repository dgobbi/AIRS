/*=========================================================================

  Module: vtkMultiLevelImageRegistration.h

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkMultiLevelImageRegistration - a multi-resolution approach.
// .SECTION Description
// This class will perform image registration by starting with low-resolution
// (blurred and subsampled) images, and progressing though several levels of
// increasing the resolution.  This accelerates the search and increases the
// probability of finding the global optimum for the registration.

#ifndef vtkMultiLevelImageRegistration_h
#define vtkMultiLevelImageRegistration_h

#include "vtkImageRegistrationBase.h"

class vtkImageRegistration;
class vtkImageResize;

class VTK_EXPORT vtkMultiLevelImageRegistration :
  public vtkImageRegistrationBase
{
public:
  vtkTypeMacro(vtkMultiLevelImageRegistration, vtkImageRegistrationBase);
  static vtkMultiLevelImageRegistration *New();
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the number of registration stages.
  // By default, each stage doubles the resolution of the previous stage.
  vtkSetMacro(NumberOfLevels, int);
  vtkGetMacro(NumberOfLevels, int);

  // Description:
  // Get the current registration level.
  vtkGetMacro(Level, int);

  // Description:
  // This is true if the current level is finished.
  bool IsLevelDone() { return this->LevelDone; }

  // Description:
  // Initialize the transform.  This will also initialize the
  // NumberOfEvaluations to zero.  If a TransformInitializer is
  // set, then only the rotation part of this matrix will be used,
  // and the initial translation will be set from the initializer.
  void Initialize(vtkMatrix4x4 *matrix);

  // Description:
  // Iterate the registration.  Returns zero if the termination condition has
  // been reached.
  int Iterate();

protected:
  vtkMultiLevelImageRegistration();
  ~vtkMultiLevelImageRegistration();

  // Description:
  // Initialize a specific level of the multi-level registration.
  // This is automatically called when necessary.
  void InitializeLevel(int level);

  int ExecuteRegistration();

  vtkImageRegistration *Helper;

  int Level;
  int NumberOfLevels;
  bool LevelDone;

private:
  // Copy constructor and assigment operator are purposely not implemented
  vtkMultiLevelImageRegistration(const vtkMultiLevelImageRegistration&);
  void operator=(const vtkMultiLevelImageRegistration&);
};

#endif /* vtkMultiLevelImageRegistration_h */
