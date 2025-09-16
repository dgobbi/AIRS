/*=========================================================================

  Program:   Atamai Classes for VTK
  Module:    vtkFrameFinder.h

Copyright (c) 2005 Atamai, Inc.
All rights reserved.

Use, modification and redistribution of the software, in source or
binary forms, are permitted provided that the following terms and
conditions are met:

1) Redistribution of the source code, in verbatim or modified
   form, must retain the above copyright notice, this license,
   the following disclaimer, and any notices that refer to this
   license and/or the following disclaimer.

2) Redistribution in binary form must include the above copyright
   notice, a copy of this license and the following disclaimer
   in the documentation or with other materials provided with the
   distribution.

3) Modified copies of the source code must be clearly marked as such,
   and must not be misrepresented as verbatim copies of the source code.

THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE SOFTWARE "AS IS"
WITHOUT EXPRESSED OR IMPLIED WARRANTY INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE.  IN NO EVENT SHALL ANY COPYRIGHT HOLDER OR OTHER PARTY WHO MAY
MODIFY AND/OR REDISTRIBUTE THE SOFTWARE UNDER THE TERMS OF THIS LICENSE
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, LOSS OF DATA OR DATA BECOMING INACCURATE
OR LOSS OF PROFIT OR BUSINESS INTERRUPTION) ARISING IN ANY WAY OUT OF
THE USE OR INABILITY TO USE THE SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES.

=========================================================================*/
// .NAME vtkFrameFinder - Locate a Leksell frame in an image.
// .SECTION Description
// This class will attempt to locate a Leksell frame within its input image.

#ifndef vtkFrameFinder_h
#define vtkFrameFinder_h

#include "vtkImageRegistrationModule.h" // For export macro
#include "vtkAlgorithm.h"

class vtkImageData;
class vtkPolyData;
class vtkMatrix4x4;

class VTKIMAGEREGISTRATION_EXPORT vtkFrameFinder : public vtkAlgorithm
{
public:
  vtkTypeMacro(vtkFrameFinder, vtkAlgorithm);
  static vtkFrameFinder *New();
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // The input to this filter should be an image of a head with a
  // stereotactic frame.  Axial images are required.  Call Update()
  // in order to compute the registration matrix.
  void SetInputData(vtkDataObject *);
  void SetInput(vtkDataObject *o) { this->SetInputData(o); }
  vtkDataObject *GetInput();

  // Description:
  // The first output of this filter is a polydata that contains all
  // of the points in the image that were identified as being part of
  // the frame fiducials.  The second output is a wireframe model of
  // the frame fiducials.  If the frame-finding algorithm succeeded,
  // then the ImageToFrameMatrix can be used to transform the found
  // points onto the frame model.
  vtkPolyData *GetOutput(int i);
  vtkPolyData *GetOutput() { return this->GetOutput(0); }

  // Description:
  // Get a boolean value that indicates whether
  bool GetSuccess() { return this->Success; }

  // Description:
  // Specify whether to use the anterior and posterior fiducials.
  // By default both of these are on, but sometimes these plates
  // are not securely fastened and therefore not reliable.  The
  // default is to use them if they are detected in the image.
  vtkGetMacro(UseAnteriorFiducial, int);
  vtkSetMacro(UseAnteriorFiducial, int);
  vtkBooleanMacro(UseAnteriorFiducial, int);
  vtkGetMacro(UsePosteriorFiducial, int);
  vtkSetMacro(UsePosteriorFiducial, int);
  vtkBooleanMacro(UsePosteriorFiducial, int);

  // Description:
  // Set the matrix that converts image data coordinates into DICOM
  // patient coordinates.  Only the first two columns of this matrix
  // are used, which must be set from the ImageOrientationPatient
  // metadata.
  void SetDICOMPatientMatrix(vtkMatrix4x4 *matrix);

  // Description:
  // After this filter has been updated, this provides the transformation
  // from image data coordinates to frame coordinates.
  vtkMatrix4x4 *GetImageToFrameMatrix() {
    return this->ImageToFrameMatrix; }

protected:
  vtkFrameFinder();
  ~vtkFrameFinder();

  // Functions overridden from Superclass
  virtual int ProcessRequest(vtkInformation *,
                             vtkInformationVector **,
                             vtkInformationVector *);
  virtual int RequestData(vtkInformation *,
			  vtkInformationVector **,
			  vtkInformationVector *);
  virtual int RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inInfo,
                                 vtkInformationVector *vtkNotUsed(outInfo));
  virtual int RequestInformation(vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inInfo,
                                 vtkInformationVector *vtkNotUsed(outInfo));
  virtual int FillInputPortInformation(int port, vtkInformation* info);
  virtual int FillOutputPortInformation(int port, vtkInformation* info);

  int FindFrame(vtkImageData *, vtkPolyData *,
                const double direction[3], vtkMatrix4x4 *matrix);

  vtkMatrix4x4 *ImageToFrameMatrix;
  vtkMatrix4x4 *DICOMPatientMatrix;

  bool Success;
  int UseAnteriorFiducial;
  int UsePosteriorFiducial;

private:
  // Copy constructor and assigment operator are purposely not implemented
  vtkFrameFinder(const vtkFrameFinder&);
  void operator=(const vtkFrameFinder&);
};

#endif //vtkFrameFinder_h
