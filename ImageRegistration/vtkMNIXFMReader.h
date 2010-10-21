/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkMNIXFMReader.h,v $

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
// .NAME vtkMNIXFMReader - A reader for MNI transformation files.
// .SECTION Description
// The MNI .xfm file format is used to store geometrical
// transformations.  Three kinds of transformations are supported by
// the file format: affine, thin-plate spline, and grid transformations.
// .SECTION See Also
// vtkMINCReader vtkMNIXFMWriter

#ifndef __vtkMNIXFMReader_h
#define __vtkMNIXFMReader_h

#include "vtkAlgorithm.h"

class vtkAbstractTransform;
class vtkDoubleArray;
class vtkCollection;

class VTK_EXPORT vtkMNIXFMReader : public vtkAlgorithm
{
public:
  vtkTypeRevisionMacro(vtkMNIXFMReader,vtkAlgorithm);

  static vtkMNIXFMReader *New();
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the file name.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Get the entension for this file format.
  virtual const char* GetFileExtensions() {
    return ".xfm"; }

  // Description:
  // Get the name of this file format.
  virtual const char* GetDescriptiveName() {
    return "MNI Transform"; }

  // Description:
  // Test whether the specified file can be read.
  virtual int CanReadFile(const char* name);

  // Description:
  // Get the number of transforms in the file.
  virtual int GetNumberOfTransforms();

  // Description:
  // Get one of the transforms listed in the file.
  virtual vtkAbstractTransform *GetNthTransform(int i);

  // Description:
  // Get the transform that results from concatenating all
  // of the transforms in the file.  This will return null
  // if you have not specified a file name.
  virtual vtkAbstractTransform *GetTransform();

  // Description:
  // Get any comments that are included in the file.
  virtual const char *GetComments();

protected:
  vtkMNIXFMReader();
  ~vtkMNIXFMReader();

  char *FileName;
  vtkAbstractTransform *Transform;
  vtkCollection *Transforms;
  int LineNumber;
  char *Comments;

  void SetTransform(vtkAbstractTransform *transform);

  int ReadLine(istream &infile, char result[256]);
  int ReadLineAfterComments(istream &infile, char result[256]);
  int SkipWhitespace(istream &infile, char linetext[256], char **cpp);
  int ParseLeftHandSide(istream &infile, char linetext[256], char **cpp,
                        char identifier[256]);
  int ParseStringValue(istream &infile, char linetext[256], char **cpp,
                       char data[256]);
  int ParseFloatValues(istream &infile, char linetext[256], char **cpp,
                       vtkDoubleArray *array);
  int ParseInvertFlagValue(istream &infile, char linetext[256], char **cpp,
                           int *invertFlag);

  int ReadLinearTransform(istream &infile, char linetext[256], char **cp);
  int ReadThinPlateSplineTransform(istream &infile, char linetext[256],
                                   char **cp);
  int ReadGridTransform(istream &infile, char linetext[256], char **cp);

  virtual int ReadNextTransform(istream &infile, char linetext[256]);

  virtual int ReadFile();

  virtual int ProcessRequest(vtkInformation* request,
                             vtkInformationVector** inInfo,
                             vtkInformationVector* outInfo);

private:
  vtkMNIXFMReader(const vtkMNIXFMReader&); // Not implemented
  void operator=(const vtkMNIXFMReader&);  // Not implemented

};

#endif
