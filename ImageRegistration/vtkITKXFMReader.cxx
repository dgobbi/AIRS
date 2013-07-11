/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkITKXFMReader.cxx

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

#include "vtkITKXFMReader.h"

#include "vtkObjectFactory.h"

#include "vtkImageData.h"
#include "vtkCollection.h"
#include "vtkTransform.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkDemandDrivenPipeline.h"
#include "vtkMath.h"
#include "vtkSmartPointer.h"

#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <string>
#include <vector>
#include <vtksys/SystemTools.hxx>

#define ITKXFM_MAXLINE 1024

//--------------------------------------------------------------------------
vtkStandardNewMacro(vtkITKXFMReader);

//-------------------------------------------------------------------------
vtkITKXFMReader::vtkITKXFMReader()
{
  this->FileName = 0;
  this->Transform = 0;
  this->Transforms = vtkCollection::New();
  this->LineNumber = 0;
  this->Comments = 0;
}

//-------------------------------------------------------------------------
vtkITKXFMReader::~vtkITKXFMReader()
{
  if (this->Transforms)
    {
    this->Transforms->Delete();
    }
  if (this->Transform)
    {
    this->Transform->Delete();
    }

  delete [] this->FileName;
  delete [] this->Comments;
}

//-------------------------------------------------------------------------
void vtkITKXFMReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "none") << "\n";
  os << indent << "Transform: " << this->Transform << "\n";
  if (this->Transform)
    {
    this->Transform->PrintSelf(os, indent.GetNextIndent());
    }
  os << indent << "NumberOfTransforms: "
     << this->Transforms->GetNumberOfItems() << "\n";
  os << indent << "Comments: "
     << (this->Comments ? this->Comments : "none") << "\n";
}

//-------------------------------------------------------------------------
int vtkITKXFMReader::CanReadFile(const char* fname)
{
  // First make sure the file exists.  This prevents an empty file
  // from being created on older compilers.
  struct stat fs;
  if(stat(fname, &fs) != 0)
    {
    return 0;
    }

  // Try to read the first line of the file.
  int status = 0;

  ifstream infile(fname);

  if (infile.good())
    {
    status = 1;
    char linetext[ITKXFM_MAXLINE];
    infile.getline(linetext, ITKXFM_MAXLINE);
    if (strncmp(linetext, "#Insight Transform File", 23) != 0)
      {
      status = 0;
      }

    infile.close();
    }

  return status;
}

//-------------------------------------------------------------------------
// Internal function to read in a line up to ITKXFM_MAXLINE characters and
// then skip to the next line in the file.
int vtkITKXFMReader::ReadLine(
  istream &infile, char result[ITKXFM_MAXLINE])
{
  this->LineNumber++;

  infile.getline(result,ITKXFM_MAXLINE);
  if (infile.fail())
    {
    if (infile.eof())
      {
      return 0;
      }
    if (infile.gcount() == ITKXFM_MAXLINE-1)
      {
      // Read ITKXFM_MAXLINE chars; ignoring the rest of the line.
      infile.clear();
      infile.ignore(VTK_INT_MAX, '\n');
      vtkWarningMacro("Overlength line (limit is " << (ITKXFM_MAXLINE-1)
                      << ") in "
                      << this->FileName << ":" << this->LineNumber);
      }
    }

  return 1;
}

//-------------------------------------------------------------------------
// Skip all blank lines or comment lines and return the first useful line,
// also set cpp to the beginning of the line for convenience.
int vtkITKXFMReader::ReadLineAfterComments(
  istream &infile, char result[ITKXFM_MAXLINE])
{
  // Skip over any comment lines or blank lines.
  // Comment lines start with '#'
  do
    {
    this->ReadLine(infile, result);
    }
  while (infile.good() && result[0] == '#');

  return 0;
}

//-------------------------------------------------------------------------
// Skip all whitespace, reading additional lines if necessary
int vtkITKXFMReader::SkipWhitespace(char **cpp)
{
  char *cp = *cpp;

  // Skip leading whitespace
  while (isspace(*cp))
    {
    cp++;
    }

  *cpp = cp;

  return (*cp != '\0');
}

//-------------------------------------------------------------------------
// Read the left hand side of a statement, including the colon
// and any whitespace following the colon.
int vtkITKXFMReader::ParseLeftHandSide(
  char **cpp, char identifier[ITKXFM_MAXLINE])
{
  int i = 0;
  char *cp = *cpp;

  // Read alphanumeric plus underscore
  if (!isdigit(*cp))
    {
    while (isalnum(*cp) || *cp == '_')
      {
      identifier[i++] = *cp++;
      }
    }
  identifier[i] = '\0';

  // Skip trailing whitespace
  while (isspace(*cp))
    {
    cp++;
    }

  // Check for equals
  this->SkipWhitespace(&cp);
  if (*cp != ':')
    {
    vtkErrorMacro("Missing \':\' " << this->FileName
                  << ":" << this->LineNumber);
    return 0;
    }
  cp++;

  // Skip ahead to the value part of the statement
  this->SkipWhitespace(&cp);

  *cpp = cp;

  return 1;
}

//-------------------------------------------------------------------------
// Read a string value.  Any leading/trailing whitespace will be ignored.
// string may not be split across multiple lines.
int vtkITKXFMReader::ParseStringValue(
  char **cpp, char data[ITKXFM_MAXLINE])
{
  char *cp = *cpp;
  this->SkipWhitespace(&cp);

  // Read until end of the line
  int i = 0;
  while (*cp)
    {
    data[i++] = *cp++;
    }

  // Remove trailing whitespace
  while (i > 0 && isspace(data[i-1]))
    {
    i--;
    }

  data[i] = '\0';

  *cpp = cp;

  return 1;
}

//-------------------------------------------------------------------------
// Read floating-point values into a vtkDoubleArray until the end of the
// line is reached.
int vtkITKXFMReader::ParseFloatValues(
  char **cpp, vtkDoubleArray *array)
{
  char *cp = *cpp;
  this->SkipWhitespace(&cp);

  // Read until end of the line
  while (*cp && !isspace(*cp))
    {
    char *tmp = cp;
    double val = strtod(cp, &cp);
    if (cp == tmp)
      {
      vtkErrorMacro("Syntax error " << this->FileName
                    << ":" << this->LineNumber);
      return 0;
      }
    array->InsertNextValue(val);
    this->SkipWhitespace(&cp);
    }

  *cpp = cp;

  return 1;
}

//-------------------------------------------------------------------------
int vtkITKXFMReader::CheckNumberOfParameters(
  vtkDoubleArray *parameters, vtkDoubleArray *fixedParameters,
  vtkIdType n, vtkIdType m, const char *classname)
{
  if (parameters->GetNumberOfTuples() != n)
    {
    vtkErrorMacro("Incorrect number of Parameters for \'"
                  << classname << "\' in "
                  << this->FileName << ":" << this->LineNumber);
    return 0;
    }
  if (fixedParameters->GetNumberOfTuples() != m)
    {
    vtkErrorMacro("Incorrect number of FixedParameters for \'"
                  << classname << "\' in "
                  << this->FileName << ":" << this->LineNumber);
    return 0;
    }

  return 1;
}

//-------------------------------------------------------------------------
void vtkITKXFMReader::BuildTransform(
  const double matparms[9], const double translation[3],
  const double center[3], vtkTransform *transform)
{
  double matrix[16] = {
    1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 1.0 };

  double mc[4] = { 0.0, 0.0, 0.0, 0.0 };
  for (int i = 0; i < 3; i++)
    {
    matrix[4*i + 0] = matparms[3*i + 0];
    matrix[4*i + 1] = matparms[3*i + 1];
    matrix[4*i + 2] = matparms[3*i + 2];
    mc[i] = center[i];
    }

  vtkMatrix4x4::MultiplyPoint(matrix, mc, mc);

  // compute the offset
  matrix[3] = translation[0] + center[0] - mc[0];
  matrix[7] = translation[1] + center[1] - mc[1];
  matrix[11] = translation[2] + center[2] - mc[2];

  transform->Concatenate(matrix);
}

//-------------------------------------------------------------------------
void vtkITKXFMReader::MatrixFromQuaternion(
  const double quat[4], double matparms[9])
{
  double mat[3][3];
  vtkMath::QuaternionToMatrix3x3(quat, mat);

  for (int i = 0; i < 3; i++)
    {
    matparms[3*i + 0] = mat[i][0];
    matparms[3*i + 1] = mat[i][1];
    matparms[3*i + 2] = mat[i][2];
    }
}

//-------------------------------------------------------------------------
void vtkITKXFMReader::MatrixFromVersor(
  const double versor[3], double matparms[9])
{
  double quat[4];
  double n = vtkMath::Norm(versor);
  double d = 1.0 - n*n;
  d = (d > 0 ? d : 0.0);
  quat[0] = versor[0];
  quat[1] = versor[1];
  quat[2] = versor[2];
  quat[3] = sqrt(d);

  vtkITKXFMReader::MatrixFromQuaternion(quat, matparms);
}

//-------------------------------------------------------------------------
void vtkITKXFMReader::MatrixFromEuler(
  const double xyz[3], double matparms[9])
{
  double cx = cos(xyz[0]);
  double sx = sin(xyz[0]);
  double cy = cos(xyz[1]);
  double sy = sin(xyz[1]);
  double cz = cos(xyz[2]);
  double sz = sin(xyz[2]);

  // rotate Y, then X, then Z
  matparms[0] = cy*cz - sx*sy*sz;
  matparms[1] = cy*sz + sx*sy*cz;
  matparms[2] = -cx*sy;

  matparms[3] = -cx*sz;
  matparms[4] = cx*cz;
  matparms[5] = sx;

  matparms[6] = sy*cz + sx*cy*sz;
  matparms[7] = sy*sz - sx*cy*cz;
  matparms[8] = cx*cy;
}

//-------------------------------------------------------------------------
int vtkITKXFMReader::ReadTransform(
  istream &infile, char linetext[ITKXFM_MAXLINE])
{
  char identifier[ITKXFM_MAXLINE];
  char classname[ITKXFM_MAXLINE];
  char *cp;

  this->ReadLineAfterComments(infile, linetext);
  cp = linetext;

  // Check for transform name
  if (this->ParseLeftHandSide(&cp, identifier) &&
      strcmp(identifier, "Transform") != 0)
    {
    vtkErrorMacro("Expected \'Transform\' in "
                  << this->FileName << ":" << this->LineNumber);
    return 0;
    }

  if (!this->ParseStringValue(&cp, classname))
    {
    vtkErrorMacro("Expected a transform type in "
                  << this->FileName << ":" << this->LineNumber);
    return 0;
    }

  // Get the transform parameters
  this->ReadLine(infile, linetext);
  cp = linetext;
  if (this->ParseLeftHandSide(&cp, identifier) &&
      strcmp(identifier, "Parameters") != 0)
    {
    vtkErrorMacro("Expected \'Parameters\' in "
                  << this->FileName << ":" << this->LineNumber);
    return 0;
    }

  vtkSmartPointer<vtkDoubleArray> parameters =
    vtkSmartPointer<vtkDoubleArray>::New();
  if (!this->ParseFloatValues(&cp, parameters))
    {
    return 0;
    }

  // Get the transform fixed parameters
  this->ReadLine(infile, linetext);
  cp = linetext;
  if (this->ParseLeftHandSide(&cp, identifier) &&
      strcmp(identifier, "FixedParameters") != 0)
    {
    vtkErrorMacro("Expected \'FixedParameters\' in "
                  << this->FileName << ":" << this->LineNumber);
    return 0;
    }

  vtkSmartPointer<vtkDoubleArray> fixedParameters =
    vtkSmartPointer<vtkDoubleArray>::New();
  if (!this->ParseFloatValues(&cp, fixedParameters))
    {
    return 0;
    }

  // Get the template parameters
  int d1 = 3;
  int d2 = 3;
  int pcount = 0;
  for (int l = 0; classname[l] != '\0'; l++)
    {
    if (classname[l] == '_')
      {
      classname[l] = '\0';
      if (pcount == 1 || pcount == 2)
        {
        int d = static_cast<int>(atol(&classname[l+1]));
        if (pcount == 1 && d != 0)
          {
          d1 = d;
          d2 = d;
          }
        else if (pcount == 2 && d != 0)
          {
          d2 = d;
          }
        }
      pcount++;
      }
    }

  vtkSmartPointer<vtkTransform> transform =
    vtkSmartPointer<vtkTransform>::New();

  double matparms[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
  double translation[3] = { 0.0, 0.0, 0.0 };
  double center[3] = { 0.0, 0.0, 0.0 };

  if (strcmp(classname, "IdentityTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, 0, 0, classname))
      {
      return 0;
      }
    }
  else if (strcmp(classname, "TranslationTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, d2, 0, classname))
      {
      return 0;
      }

    int rows = (d2 < 3 ? d2 : 3);
    for (int i = 0; i < rows; i++)
      {
      translation[i] = parameters->GetValue(i);
      }
    }
  else if (strcmp(classname, "ScaleTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, d2, 0, classname))
      {
      return 0;
      }

    int rows = (d2 < 3 ? d2 : 3);
    for (int i = 0; i < rows; i++)
      {
      matparms[3*i + i] = parameters->GetValue(i);
      }
    }
  else if (strcmp(classname, "ScaleLogarithmicTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, d2, 0, classname))
      {
      return 0;
      }

    int rows = (d2 < 3 ? d2 : 3);
    for (int i = 0; i < rows; i++)
      {
      matparms[3*i + i] = exp(parameters->GetValue(i));
      }
    }
  else if (strcmp(classname, "MatrixOffsetTransformBase") == 0 ||
           strcmp(classname, "AffineTransform") == 0 ||
           strcmp(classname, "Rigid3DTransform") == 0 ||
           strcmp(classname, "v3Rigid3DTransform") == 0 ||
           strcmp(classname, "ScalableAffineTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, (d1+1)*d2, d2, classname))
      {
      return 0;
      }

    int rows = (d2 < 3 ? d2 : 3);
    int cols = (d1 < 3 ? d1 : 3);
    for (int i = 0; i < rows; i++)
      {
      for (int j = 0; j < cols; j++)
        {
        matparms[3*i + j] = parameters->GetValue(d1*i + j);
        }
      translation[i] = parameters->GetValue(d1*d2 + i);
      }

    for (int j = 0; j < cols; j++)
      {
      center[j] = fixedParameters->GetValue(j);
      }
    }
  else if (strcmp(classname, "CenteredAffineTransform") == 0 ||
           strcmp(classname, "FixedCenterOfRotationAffineTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, (d1+1)*d2 + d1, d2, classname))
      {
      return 0;
      }

    int rows = (d2 < 3 ? d2 : 3);
    int cols = (d1 < 3 ? d1 : 3);
    for (int i = 0; i < rows; i++)
      {
      for (int j = 0; j < cols; j++)
        {
        matparms[3*i + j] = parameters->GetValue(d1*i + j);
        }
      translation[i] = parameters->GetValue(d1*(d2 + 1) + i);
      }

    for (int j = 0; j < cols; j++)
      {
      center[j] = parameters->GetValue(d1*d2 + j);
      }
    }
  else if (strcmp(classname, "Euler3DTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, 6, 3, classname))
      {
      return 0;
      }

    double xyz[3];
    for (int i = 0; i < 3; i++)
      {
      xyz[i] = parameters->GetValue(i);
      translation[i] = parameters->GetValue(3 + i);
      center[i] = fixedParameters->GetValue(i);
      }
    vtkITKXFMReader::MatrixFromEuler(xyz, matparms);
    }
  else if (strcmp(classname, "Rigid2DTransform") == 0 ||
           strcmp(classname, "Euler2DTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, 3, 2, classname))
      {
      return 0;
      }

    double angle = parameters->GetValue(0);
    translation[0] = parameters->GetValue(1);
    translation[1] = parameters->GetValue(2);
    center[0] = fixedParameters->GetValue(0);
    center[1] = fixedParameters->GetValue(1);
    matparms[0] = cos(angle);
    matparms[1] = -sin(angle);
    matparms[3] = sin(angle);
    matparms[4] = cos(angle);
    }
  else if (strcmp(classname, "Similarity2DTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, 4, 2, classname))
      {
      return 0;
      }

    double scale = parameters->GetValue(0);
    double angle = parameters->GetValue(1);
    translation[0] = parameters->GetValue(2);
    translation[1] = parameters->GetValue(3);
    center[0] = fixedParameters->GetValue(0);
    center[1] = fixedParameters->GetValue(1);
    matparms[0] = cos(angle)*scale;
    matparms[1] = -sin(angle)*scale;
    matparms[3] = sin(angle)*scale;
    matparms[4] = cos(angle)*scale;
    }
  else if (strcmp(classname, "CenteredRigid2DTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, 5, 2, classname))
      {
      return 0;
      }

    double angle = parameters->GetValue(0);
    center[0] = fixedParameters->GetValue(1);
    center[1] = fixedParameters->GetValue(2);
    translation[0] = parameters->GetValue(3);
    translation[1] = parameters->GetValue(4);
    matparms[0] = cos(angle);
    matparms[1] = -sin(angle);
    matparms[3] = sin(angle);
    matparms[4] = cos(angle);
    }
  else if (strcmp(classname, "CenteredSimilarity2DTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, 6, 2, classname))
      {
      return 0;
      }

    double scale = parameters->GetValue(0);
    double angle = parameters->GetValue(1);
    center[0] = parameters->GetValue(2);
    center[1] = parameters->GetValue(3);
    translation[0] = parameters->GetValue(4);
    translation[1] = parameters->GetValue(5);
    matparms[0] = cos(angle)*scale;
    matparms[1] = -sin(angle)*scale;
    matparms[3] = sin(angle)*scale;
    matparms[4] = cos(angle)*scale;
    }
  else if (strcmp(classname, "CenteredEuler3DTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, 9, 3, classname))
      {
      return 0;
      }

    double xyz[3];
    for (int i = 0; i < 3; i++)
      {
      xyz[i] = parameters->GetValue(i);
      center[i] = parameters->GetValue(3 + i);
      translation[i] = parameters->GetValue(6 + i);
      }
    vtkITKXFMReader::MatrixFromEuler(xyz, matparms);
    }
  else if (strcmp(classname, "VersorTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, 3, 3, classname))
      {
      return 0;
      }

    double versor[3];
    for (int i = 0; i < 3; i++)
      {
      versor[i] = parameters->GetValue(i);
      center[i] = fixedParameters->GetValue(i);
      }
    vtkITKXFMReader::MatrixFromVersor(versor, matparms);
    }
  else if (strcmp(classname, "VersorRigid3DTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, 6, 3, classname))
      {
      return 0;
      }

    double versor[3];
    for (int i = 0; i < 3; i++)
      {
      versor[i] = parameters->GetValue(i);
      translation[i] = parameters->GetValue(3 + i);
      center[i] = fixedParameters->GetValue(i);
      }
    vtkITKXFMReader::MatrixFromVersor(versor, matparms);
    }
  else if (strcmp(classname, "Similarity3DTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, 7, 3, classname))
      {
      return 0;
      }

    double versor[3];
    for (int i = 0; i < 3; i++)
      {
      versor[i] = parameters->GetValue(i);
      translation[i] = parameters->GetValue(3 + i);
      center[i] = fixedParameters->GetValue(i);
      }
    vtkITKXFMReader::MatrixFromVersor(versor, matparms);
    double scale = parameters->GetValue(6);
    for (int j = 0; j < 9; j++)
      {
      matparms[j] *= scale;
      }
    }
  else if (strcmp(classname, "ScaleVersor3DTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, 9, 3, classname))
      {
      return 0;
      }

    double versor[3];
    double scale[3];
    for (int i = 0; i < 3; i++)
      {
      versor[i] = parameters->GetValue(i);
      translation[i] = parameters->GetValue(3 + i);
      scale[i] = parameters->GetValue(6 + i);
      center[i] = fixedParameters->GetValue(i);
      }
    vtkITKXFMReader::MatrixFromVersor(versor, matparms);
    // just duplicating what ITK does...
    matparms[0] += scale[0] - 1.0;
    matparms[4] += scale[1] - 1.0;
    matparms[8] += scale[2] - 1.0;
    }
  else if (strcmp(classname, "ScaleSkewVersor3DTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, 15, 3, classname))
      {
      return 0;
      }

    double versor[3];
    double scale[3];
    double skew[6];
    for (int i = 0; i < 3; i++)
      {
      versor[i] = parameters->GetValue(i);
      translation[i] = parameters->GetValue(3 + i);
      scale[i] = parameters->GetValue(6 + i);
      skew[i] = parameters->GetValue(9 + i);
      skew[i+3] = parameters->GetValue(12 + i);
      center[i] = fixedParameters->GetValue(i);
      }
    vtkITKXFMReader::MatrixFromVersor(versor, matparms);
    matparms[0] += scale[0] - 1.0;
    matparms[1] += skew[0];
    matparms[2] += skew[1];
    matparms[3] += skew[2];
    matparms[4] += scale[1] - 1.0;
    matparms[5] += skew[3];
    matparms[6] += skew[4];
    matparms[7] += skew[5];
    matparms[8] += scale[2] - 1.0;
    }
  else if (strcmp(classname, "QuaternionRigidTransform") == 0)
    {
    if (!this->CheckNumberOfParameters(
          parameters, fixedParameters, 7, 3, classname))
      {
      return 0;
      }

    double quat[4];
    quat[0] = parameters->GetValue(3);

    for (int i = 0; i < 3; i++)
      {
      quat[i+1] = parameters->GetValue(i);
      translation[i] = parameters->GetValue(4 + i);
      center[i] = fixedParameters->GetValue(i);
      }
    vtkITKXFMReader::MatrixFromQuaternion(quat, matparms);
    }

  else
    {
    vtkErrorMacro("Unrecognized transform type \'"
                  << classname << "\' in ");
    return 0;
    }

  vtkITKXFMReader::BuildTransform(
    matparms, translation, center, transform);

  this->Transforms->AddItem(transform);

  return 1;
}

//-------------------------------------------------------------------------
int vtkITKXFMReader::ReadFile()
{
  this->Transforms->RemoveAllItems();
  this->SetTransform(0);

  // Check that the file name has been set.
  if (!this->FileName || this->FileName[0] == '\0')
    {
    vtkErrorMacro("ReadFile: No file name has been set");
    return 0;
    }

  // Make sure that the file exists.
  struct stat fs;
  if(stat(this->FileName, &fs) != 0)
    {
    vtkErrorMacro("ReadFile: Can't open file " << this->FileName);
    return 0;
    }

  // Make sure that the file is readable.
  ifstream infile(this->FileName);

  if (infile.fail())
    {
    vtkErrorMacro("ReadFile: Can't read the file " << this->FileName);
    return 0;
    }

  // Read the first line
  char linetext[ITKXFM_MAXLINE];
  this->LineNumber = 0;
  this->ReadLine(infile, linetext);

  if (strncmp(linetext, "#Insight Transform File", 23) != 0)
    {
    vtkErrorMacro("ReadFile: File is not an ITK transform file: " << this->FileName);
    infile.close();
    return 0;
    }
  else
    {
    size_t n = strlen(linetext);
    while (n > 1 && isspace(linetext[n-1]))
      {
      --n;
      }
    delete [] this->Comments;
    this->Comments = new char[n + 1];
    strncpy(this->Comments, linetext, n);
    this->Comments[n] = '\0';
    }

  // Read the transforms
  while (infile.good())
    {
    if (this->ReadTransform(infile, linetext) == 0)
      {
      this->Transforms->RemoveAllItems();
      infile.close();
      return 0;
      }
    this->ReadLine(infile, linetext);
    }

  // Close the file
  infile.close();

  // Create the output transform.
  int n = this->Transforms->GetNumberOfItems();
  if (n == 1)
    {
    this->SetTransform(
      (vtkTransform *)this->Transforms->GetItemAsObject(0));
    }
  else
    {
    // Determine whether the full transform is linear
    int linear = 1;
    int i = 0;
    for (i = 0; i < n; i++)
      {
      if (!this->Transforms->GetItemAsObject(i)->IsA("vtkLinearTransform"))
        {
        linear = 0;
        break;
        }
      }

    // If linear, use vtkTransform to concatenate,
    // else use vtkGeneralTransform.
    if (linear)
      {
      vtkTransform *transform = vtkTransform::New();
      transform->PostMultiply();
      for (i = 0; i < n; i++)
        {
        vtkLinearTransform *linearTransform =
          (vtkLinearTransform *)this->Transforms->GetItemAsObject(i);
        transform->Concatenate(linearTransform->GetMatrix());
        }
      this->SetTransform(transform);
      transform->Delete();
      }
    else
      {
      vtkTransform *transform = vtkTransform::New();
      transform->PostMultiply();
      for (i = 0; i < n; i++)
        {
        vtkTransform *aTransform =
          (vtkTransform *)this->Transforms->GetItemAsObject(i);
        transform->Concatenate(aTransform->GetMatrix());
        }
      this->SetTransform(transform);
      transform->Delete();
      }
    }

  return 1;
}

//-------------------------------------------------------------------------
int vtkITKXFMReader::ProcessRequest(vtkInformation *request,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector)
{
  if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
    {
    return this->ReadFile();
    }

  return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

//-------------------------------------------------------------------------
void vtkITKXFMReader::SetTransform(vtkTransform *transform)
{
  if (this->Transform != transform)
    {
    if (this->Transform)
      {
      this->Transform->Delete();
      }
    if (transform)
      {
      transform->Register(this);
      }
    this->Transform = transform;
    }
}

//-------------------------------------------------------------------------
vtkTransform *vtkITKXFMReader::GetTransform()
{
  this->Update();

  return this->Transform;
}

//-------------------------------------------------------------------------
int vtkITKXFMReader::GetNumberOfTransforms()
{
  this->Update();

  return this->Transforms->GetNumberOfItems();
}

//-------------------------------------------------------------------------
vtkTransform *vtkITKXFMReader::GetNthTransform(int i)
{
  this->Update();

  if (i < 0 || i >= this->Transforms->GetNumberOfItems())
    {
    return 0;
    }

  return (vtkTransform *)this->Transforms->GetItemAsObject(i);
}

//-------------------------------------------------------------------------
const char *vtkITKXFMReader::GetComments()
{
  this->Update();

  if (this->Comments == 0)
    {
    return "";
    }

  return this->Comments;
}
