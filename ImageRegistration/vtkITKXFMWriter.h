/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkITKXFMWriter.h

Copyright (c) 2006 Atamai, Inc.

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
// .NAME vtkITKXFMWriter - A writer for ITK transformation files.
// .SECTION Description
// ITK uses a text file to save its transformations.  The file contains
// the name of the transform class and the parameters.  Not all ITK
// transform types are supported by this reader, only linear transforms
// that can be represented by a 4x4 matrix.
// .SECTION See Also
// vtkITKTransformReader vtkMNITransformReader
// .SECTION Thanks
// Thanks to David Gobbi for writing this class.

#ifndef __vtkITKXFMWriter_h
#define __vtkITKXFMWriter_h

#include "vtkImageRegistrationModule.h" // For export macro
#include "vtkAlgorithm.h"

class vtkAbstractTransform;
class vtkHomogeneousTransform;
class vtkThinPlateSplineTransform;
class vtkGridTransform;
class vtkCollection;

class VTKIMAGEREGISTRATION_EXPORT vtkITKXFMWriter : public vtkAlgorithm
{
public:
  vtkTypeMacro(vtkITKXFMWriter,vtkAlgorithm);

  static vtkITKXFMWriter *New();
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the file name.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Get the entension for this file format.
  virtual const char* GetFileExtensions() {
    return ".txt" ".tfm" ".mat"; }

  // Description:
  // Get the name of this file format.
  virtual const char* GetDescriptiveName() {
    return "ITK Transform"; }

  // Description:
  // Set the transform.
  virtual void SetTransform(vtkAbstractTransform *transform);
  virtual vtkAbstractTransform *GetTransform() {
    return this->Transform; };

  // Description:
  // Set the transform center.
  virtual void SetTransformCenter(double x, double y, double z);
  virtual void SetTransformCenter(const double c[3]) {
    this->SetTransformCenter(c[0], c[1], c[2]); }
  vtkGetVector3Macro(TransformCenter, double);

  // Description:
  // Add another transform to the file.  The next time that
  // SetTransform is called, all added transforms will be
  // removed.
  virtual void AddTransform(vtkAbstractTransform *transform);

  // Description:
  // Get the number of transforms that will be written.
  virtual int GetNumberOfTransforms();

  // Description:
  // Write the file.
  virtual void Write();

protected:
  vtkITKXFMWriter();
  ~vtkITKXFMWriter();

  char *FileName;
  vtkAbstractTransform *Transform;
  vtkCollection *Transforms;
  double TransformCenter[3];

  int WriteLinearTransform(
    ostream &outfile, vtkHomogeneousTransform *transform);

  virtual int WriteTransform(ostream &outfile,
                             vtkAbstractTransform *transform);

  virtual int WriteFile();

  static bool IsMatFile(const char *fname);

  virtual int ProcessRequest(vtkInformation* request,
                             vtkInformationVector** inInfo,
                             vtkInformationVector* outInfo);

private:
  vtkITKXFMWriter(const vtkITKXFMWriter&); // Not implemented
  void operator=(const vtkITKXFMWriter&);  // Not implemented

};

#endif
