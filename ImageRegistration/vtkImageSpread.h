/*=========================================================================

  Module: vtkImageSpread.h

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageSpread - Spread an image to fill outside areas
// .SECTION Description
// Given a stencil that marks "valid" regions of the image, the "invalid"
// voxels are filled with the average of any valid neighbors.  Neighbors
// are determined by face connectivity for voxels, and edge connectivity
// for pixels.  This spreading process is repeated for the specified number
// of iterations, or until all of the invalid voxels have been replaced.

#ifndef vtkImageSpread_h
#define vtkImageSpread_h

#include "vtkImageAlgorithm.h"

class vtkImageStencilData;

class VTK_EXPORT vtkImageSpread : public vtkImageAlgorithm
{
public:
  static vtkImageSpread *New();
  vtkTypeMacro(vtkImageSpread, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // The input for a stencil (input port 2).
  // Everything inside the stencil will remain unchanged, and everything
  // outside the stencil will be replaced.
  void SetStencilConnection(vtkAlgorithmOutput *port);
  vtkAlgorithmOutput *GetStencilConnection();
  void SetStencilData(vtkImageStencilData *data);

  // Description:
  // Set the number of iterations.
  vtkSetMacro(NumberOfIterations, int);
  vtkGetMacro(NumberOfIterations, int);

protected:
  vtkImageSpread();
  ~vtkImageSpread();

  int FillInputPortInformation(int port, vtkInformation *info);
  int RequestUpdateExtent(
    vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestData(
    vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  int NumberOfIterations;

private:
  vtkImageSpread(const vtkImageSpread&);  // Not implemented.
  void operator=(const vtkImageSpread&);  // Not implemented.
};

#endif
