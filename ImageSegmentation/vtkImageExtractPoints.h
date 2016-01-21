/*=========================================================================

  Copyright (c) 2015,2016 David Gobbi
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

  * Neither the name of David Gobbi nor the names of any contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.

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
// .NAME vtkImageExtractPoints - Extract all voxels within stencil
// .SECTION Description
// vtkImageExtractPoints takes an input image and a stencil, and outputs
// a flat array that contains all of the voxels withing the stencil.

#ifndef __vtkImageExtractPoints_h
#define __vtkImageExtractPoints_h

#include "vtkThreadedImageAlgorithm.h"

class vtkImageStencilData;
class vtkPoints;
class vtkIdTypeArray;
class vtkDoubleArray;

class VTK_EXPORT vtkImageExtractPoints : public vtkImageAlgorithm
{
public:
  static vtkImageExtractPoints *New();
  vtkTypeMacro(vtkImageExtractPoints,vtkImageAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the scalar component to extract from multi-component data.
  // The default value is -1, which extracts all components into a
  // multi-component array.
  vtkSetMacro(ActiveComponent, int);
  vtkGetMacro(ActiveComponent, int);

  // Description:
  // Use a stencil to compute the histogram for just a part of the image.
  void SetStencilConnection(vtkAlgorithmOutput *port);
  vtkAlgorithmOutput *GetStencilConnection();
  void SetStencilData(vtkImageStencilData *stencil);

  // Description:
  // Get all the points within the stencil.
  vtkPoints *GetPoints();

  // Description:
  // Get an array that contains all the point Ids within the stencil.
  vtkIdTypeArray *GetPointIds();

  // Description:
  // Get an array that contains all the scalars within the stencil.
  vtkDataArray *GetScalars();

protected:
  vtkImageExtractPoints();
  ~vtkImageExtractPoints();

  virtual int RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inInfo,
                                 vtkInformationVector *vtkNotUsed(outInfo));
  virtual int RequestInformation(vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inInfo,
                                 vtkInformationVector *vtkNotUsed(outInfo));
  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *);

  virtual int FillInputPortInformation(int port, vtkInformation *info);
  virtual int FillOutputPortInformation(int port, vtkInformation *info);

  int ActiveComponent;
  vtkPoints *Points;
  vtkIdTypeArray *PointIds;
  vtkDoubleArray *Scalars;

private:
  vtkImageExtractPoints(const vtkImageExtractPoints&);  // Not implemented.
  void operator=(const vtkImageExtractPoints&);  // Not implemented.
};

#endif
