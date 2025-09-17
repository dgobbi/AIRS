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
// .NAME vtkImageExtractPoints - Extract all image voxels as points.
// .SECTION Description
// This filter takes an input image and an optional stencil, and creates
// a vtkPolyData that contains the points and the point attributes but no
// cells.  If a stencil is provided, only the points inside the stencil
// are included.
// .SECTION Thanks
// Thanks to David Gobbi, Calgary Image Processing and Analysis Centre,
// University of Calgary, for providing this class.

#ifndef vtkImageExtractPoints_h
#define vtkImageExtractPoints_h

#include "vtkImageSegmentationModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"

class vtkImageStencilData;

class VTKIMAGESEGMENTATION_EXPORT vtkImageExtractPoints :
  public vtkPolyDataAlgorithm
{
public:
  static vtkImageExtractPoints *New();
  vtkTypeMacro(vtkImageExtractPoints,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Description:
  // Only extract the points that lie within the stencil.
  void SetStencilConnection(vtkAlgorithmOutput *port);
  vtkAlgorithmOutput *GetStencilConnection();
  void SetStencilData(vtkImageStencilData *stencil);

  // Description:
  // Set the desired precision for the output points.
  // See vtkAlgorithm::DesiredOutputPrecision for the available choices.
  // The default is double precision.
  vtkSetMacro(OutputPointsPrecision, int);
  vtkGetMacro(OutputPointsPrecision, int);

protected:
  vtkImageExtractPoints();
  ~vtkImageExtractPoints();

  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inInfo,
                         vtkInformationVector *outInfo) override;

  int RequestUpdateExtent(vtkInformation *request,
                          vtkInformationVector **inInfo,
                          vtkInformationVector *outInfo) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inInfo,
                  vtkInformationVector *outInfo) override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int OutputPointsPrecision;

private:
  vtkImageExtractPoints(const vtkImageExtractPoints&);  // Not implemented.
  void operator=(const vtkImageExtractPoints&);  // Not implemented.
};

#endif
