/*=========================================================================

  Program:   Atamai Image Registration and Segmentation
  Module:    ComputeStrain.cxx

  Copyright (c) 2013 David Gobbi
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

  * Neither the name of David Gobbi, nor the names of any authors nor
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.

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

#include "vtkNIFTIReader.h"
#include "vtkNIFTIWriter.h"
#include "vtkITKXFMReader.h"
#include "vtkTransformToStrain.h"
//#include "vtkImageGaussianInterpolator.h"

#include <vtkMatrix4x4.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkStringArray.h>
#include <vtkIntArray.h>
#include <vtkTransform.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkMNITransformReader.h>
#include <vtkGeneralTransform.h>
#include <vtkBSplineTransform.h>
#include <vtkImageBSplineCoefficients.h>
#include <vtkImageResize.h>
#include <vtkErrorCode.h>
#include <vtkSmartPointer.h>

#include <vtksys/SystemTools.hxx>

#include <string>
#include <vector>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

// Print the options
void strain_usage(FILE *file, const char *command_name)
{
  const char *cp = command_name + strlen(command_name);
  while (cp != command_name && cp[-1] != '\\' && cp[-1] != '/') { --cp; }

  fprintf(file,
    "usage: %s -o Output.nii [--deformation] -R Like.nii Warp.nii Affine.txt\n", cp);
  fprintf(file,
    "usage: %s -o Output.nii [--deformation] -R Like.nii -i Affine.txt InverseWarp.nii\n", cp);
  fprintf(file,
    "options:\n"
    "  -o <output.nii[.gz]>    The output file.\n"
    "  -R <file.nii[.gz]>      A file that has the size and spacing that are desired for the output.\n"
    "  --deformation           Output the deformation gradient tensor, rather than the strain tensor.\n"
    "  --help                  Print minimal documentation.\n"
  );
}

// Print the help
void strain_help(FILE *file, const char *command_name)
{
  strain_usage(file, command_name);

  fprintf(file,
    "\n");

  fprintf(file,
    "This program will compute strain tensors from a warp grid.\n"
    "\n");
  fprintf(file,
    "As input, it expects the image for which to produce the strain\n"
    "map (this image must be preceded by \"-R\"), as well as one or\n"
    "more transformations, each of which may be preceded by \"-i\" to\n"
    "indicate that the transformation is to be inverted.\n"
    "\n");
  fprintf(file,
    "As output, it produces a file that contains the strain map.  The\n"
    "output file must be preceded by \"-o\"\n."
    "\n");
}

// Print error
void strain_check_error(vtkObject *o)
{
  vtkNIFTIReader *reader = vtkNIFTIReader::SafeDownCast(o);
  vtkNIFTIWriter *writer = vtkNIFTIWriter::SafeDownCast(o);
  vtkMNITransformReader *xfmreader = vtkMNITransformReader::SafeDownCast(o);
  vtkITKXFMReader *itkreader = vtkITKXFMReader::SafeDownCast(o);
  const char *filename = 0;
  unsigned long errorcode = 0;
  if (writer)
    {
    filename = writer->GetFileName();
    errorcode = writer->GetErrorCode();
    }
  else if (reader)
    {
    filename = reader->GetInternalFileName();
    errorcode = reader->GetErrorCode();
    }
  else if (xfmreader)
    {
    filename = xfmreader->GetFileName();
    errorcode = xfmreader->GetErrorCode();
    }
  else if (itkreader)
    {
    filename = itkreader->GetFileName();
    errorcode = itkreader->GetErrorCode();
    }
  if (!filename)
    {
    filename = "";
    }

  switch(errorcode)
    {
    case vtkErrorCode::NoError:
      return;
    case vtkErrorCode::FileNotFoundError:
      fprintf(stderr, "File not found: %s\n", filename);
      break;
    case vtkErrorCode::CannotOpenFileError:
      fprintf(stderr, "Cannot open file: %s\n", filename);
      break;
    case vtkErrorCode::UnrecognizedFileTypeError:
      fprintf(stderr, "Unrecognized file type: %s\n", filename);
      break;
    case vtkErrorCode::PrematureEndOfFileError:
      fprintf(stderr, "File is truncated: %s\n", filename);
      break;
    case vtkErrorCode::FileFormatError:
      fprintf(stderr, "Bad file: %s\n", filename);
      break;
    case vtkErrorCode::NoFileNameError:
      fprintf(stderr, "Output filename could not be used: %s\n", filename);
      break;
    case vtkErrorCode::OutOfDiskSpaceError:
      fprintf(stderr, "Out of disk space while writing file: %s\n", filename);
      break;
    default:
      fprintf(stderr, "An unknown error occurred.\n");
      break;
    }

  exit(1);
}

int strain_read_transform(
  vtkGeneralTransform *transform, const char *file, bool invert,
  const double outputSpacing[3])
{
  int n = strlen(file);
  while (n) { if (file[--n] == '.') { break; } }
  if (strcmp(file + n, ".gz") == 0)
    {
    while (n) { if (file[--n] == '.') { break; } }
    }
  const char *ext = file + n;

  vtkSmartPointer<vtkAbstractTransform> t;

  if (strcmp(ext, ".xfm") == 0)
    {
    vtkSmartPointer<vtkMNITransformReader> reader =
      vtkSmartPointer<vtkMNITransformReader>::New();
    reader->SetFileName(file);
    strain_check_error(reader);
    t = reader->GetTransform();
    }
  else if (strcmp(ext, ".txt") == 0 ||
           strcmp(ext, ".tfm") == 0)
    {
    // convert ITK transforms from LPS to RAS coordinates
    // for use with NIFTI files
    static const double lps[16] = {
      -1.0, 0.0, 0.0, 0.0,
      0.0, -1.0, 0.0, 0.0,
      0.0, 0.0, 1.0, 0.0,
      0.0, 0.0, 0.0, 1.0 };
    vtkSmartPointer<vtkITKXFMReader> reader =
      vtkSmartPointer<vtkITKXFMReader>::New();
    reader->SetFileName(file);
    strain_check_error(reader);
    t = reader->GetTransform();
    vtkSmartPointer<vtkTransform> lt =
      vtkSmartPointer<vtkTransform>::New();
    lt->Concatenate(vtkLinearTransform::SafeDownCast(t)->GetMatrix());
    lt->PreMultiply();
    lt->Concatenate(lps);
    lt->PostMultiply();
    lt->Concatenate(lps);
    t = lt;
    }
  else if (strcmp(ext, ".nii") == 0 ||
           strcmp(ext, ".nii.gz") == 0)
    {
    vtkSmartPointer<vtkNIFTIReader> reader =
      vtkSmartPointer<vtkNIFTIReader>::New();
    reader->SetFileName(file);
    reader->Update();
    strain_check_error(reader);

    // the gaussian standard deviation will be 0.399 times the
    // ratio of the output image spacing to the displacement grid
    // spacing, in order to avoid aliasing of the output data
    double spacing[3], blurFactors[3];
    reader->GetOutput()->GetSpacing(spacing);
    blurFactors[0] = outputSpacing[0]/spacing[0];
    blurFactors[1] = outputSpacing[1]/spacing[1];
    blurFactors[2] = outputSpacing[2]/spacing[2];

    /*
    vtkSmartPointer<vtkImageGaussianInterpolator> interp =
      vtkSmartPointer<vtkImageGaussianInterpolator>::New();
    interp->SetBlurFactors(blurFactors);

    // smooth the grid with a Gaussian
    vtkSmartPointer<vtkImageResize> smooth =
      vtkSmartPointer<vtkImageResize>::New();
    smooth->SetInputConnection(reader->GetOutputPort());
    smooth->SetInterpolator(interp);
    smooth->Update();
    */
    vtkSmartPointer<vtkImageGaussianSmooth> smooth =
      vtkSmartPointer<vtkImageGaussianSmooth>::New();
    smooth->SetInputConnection(reader->GetOutputPort());
    smooth->SetRadiusFactors(4.5, 4.5, 4.5);
    smooth->SetStandardDeviation(
      0.399*blurFactors[0], 0.399*blurFactors[1], 0.399*blurFactors[2]);
    smooth->Update();

    // compute the b-spline coefficients
    vtkSmartPointer<vtkImageBSplineCoefficients> bsplineCoeffs =
      vtkSmartPointer<vtkImageBSplineCoefficients>::New();
    bsplineCoeffs->SetInputConnection(smooth->GetOutputPort());
    bsplineCoeffs->Update();

    // break the pipeline connection
    vtkSmartPointer<vtkImageData> image =
      vtkSmartPointer<vtkImageData>::New();
    image->CopyStructure(smooth->GetOutput());
    image->GetPointData()->PassData(smooth->GetOutput()->GetPointData());

    // reverse x and y vector components, because ITK uses LPS
    // coordinates instead of RAS like NIFTI does
    vtkDataArray *scalars = image->GetPointData()->GetScalars();
    vtkIdType m = scalars->GetNumberOfTuples();
    for (vtkIdType j = 0; j < m; j++)
      {
      double v[3];
      scalars->GetTuple(j, v);
      v[0] = -v[0];
      v[1] = -v[1];
      scalars->SetTuple(j, v);
      }

    // create a b-spline for the vectors, nice for derivative computation
    vtkSmartPointer<vtkBSplineTransform> gt =
      vtkSmartPointer<vtkBSplineTransform>::New();
    gt->SetBorderModeToZero();
    gt->SetCoefficients(image);
    t = gt;
    }
  else
    {
    fprintf(stderr, "Unrecognized transform file type \"%s\".\n", ext);
    return 0;
    }

  vtkLinearTransform *lt = vtkLinearTransform::SafeDownCast(t);
  if (lt)
    {
    vtkSmartPointer<vtkMatrix4x4> matrix =
      vtkSmartPointer<vtkMatrix4x4>::New();
    matrix->DeepCopy(lt->GetMatrix());
    if (invert)
      {
      matrix->Invert();
      }
    transform->Concatenate(matrix);
    }
  else
    {
    if (invert)
      {
      transform->Concatenate(t->GetInverse());
      }
    else
      {
      transform->Concatenate(t);
      }
    }

  return 1;
}

int main(int argc, char *argv[])
{
  // this is used to join all the transforms
  vtkSmartPointer<vtkGeneralTransform> transform =
    vtkSmartPointer<vtkGeneralTransform>::New();
  transform->PostMultiply();

  vtkSmartPointer<vtkStringArray> transformFiles =
    vtkSmartPointer<vtkStringArray>::New();
  vtkSmartPointer<vtkIntArray> transformInvert =
    vtkSmartPointer<vtkIntArray>::New();

  // the output and target file names
  const char *outfile = 0;
  const char *deffile = 0;
  const char *targetfile = 0;
  int deformation = 0;

  // read the options from the command line
  int argi = 1;
  while (argi < argc)
    {
    const char *arg = argv[argi++];
    if (strcmp(arg, "--help") == 0)
      {
      strain_help(stdout, argv[0]);
      exit(0);
      }
    else if (strcmp(arg, "--deformation") == 0)
      {
      deformation = 1;
      }
    else if (strcmp(arg, "-o") == 0)
      {
      if (argi >= argc)
        {
        fprintf(stderr, "\nA file must follow the \'-o\' flag\n\n");
        strain_usage(stderr, argv[0]);
        exit(1);
        }
      arg = argv[argi++];
      outfile = arg;
      }
    else if (strcmp(arg, "-d") == 0)
      {
      if (argi >= argc)
        {
        fprintf(stderr, "\nA file must follow the \'-d\' flag\n\n");
        strain_usage(stderr, argv[0]);
        exit(1);
        }
      arg = argv[argi++];
      deffile = arg;
      }
    else if (strcmp(arg, "-R") == 0)
      {
      if (argi >= argc)
        {
        fprintf(stderr, "\nA file must follow the \'-R\' flag\n\n");
        strain_usage(stderr, argv[0]);
        exit(1);
        }
      arg = argv[argi++];
      targetfile = arg;
      }
    else if (strcmp(arg, "-i") == 0)
      {
      if (argi >= argc)
        {
        fprintf(stderr, "\nA file must follow the \'-i\' flag\n\n");
        strain_usage(stderr, argv[0]);
        exit(1);
        }
      arg = argv[argi++];
      transformFiles->InsertNextValue(arg);
      transformInvert->InsertNextValue(1);
      }
    else if (arg[0] != '-')
      {
      transformFiles->InsertNextValue(arg);
      transformInvert->InsertNextValue(0);
      }
    else
      {
      fprintf(stderr, "\nUnrecognized option %s.\n\n", arg);
      strain_usage(stderr, argv[0]);
      exit(1);
      }
    }

  if (!outfile)
    {
    fprintf(stderr, "\nAn output file must be specified with \"-o\"\n");
    exit(1);
    }

  if (!targetfile)
    {
    fprintf(stderr, "\nA target file must be specified with \"-R\"\n");
    exit(1);
    }

  vtkSmartPointer<vtkNIFTIReader> reader =
    vtkSmartPointer<vtkNIFTIReader>::New();
  reader->SetFileName(targetfile);
  reader->Update();
  strain_check_error(reader);

  double spacing[3], origin[3];
  int extent[6];
  vtkImageData *image = reader->GetOutput();
  image->GetSpacing(spacing);
  image->GetOrigin(origin);
  image->GetExtent(extent);

  vtkIdType n = transformFiles->GetNumberOfValues();
  for (vtkIdType i = 0; i < n; i++)
    {
    std::string tfile = transformFiles->GetValue(i);
    int invert = transformInvert->GetValue(i);
    if (!strain_read_transform(transform, tfile.c_str(), invert, spacing))
      {
      exit(1);
      }
    }

  vtkSmartPointer<vtkTransformToStrain> computeStrain =
    vtkSmartPointer<vtkTransformToStrain>::New();
  if (deformation)
    {
    computeStrain->SetOutputValueToDeformationGradient();
    }
  computeStrain->SetOutputScalarTypeToFloat();
  computeStrain->SetInput(transform);
  computeStrain->SetOutputSpacing(spacing);
  computeStrain->SetOutputOrigin(origin);
  computeStrain->SetOutputExtent(extent);
  computeStrain->Update();

  vtkSmartPointer<vtkNIFTIWriter> writer =
    vtkSmartPointer<vtkNIFTIWriter>::New();
  writer->SetInputConnection(computeStrain->GetOutputPort());
  writer->SetFileName(outfile);
  writer->Write();
  strain_check_error(writer);
}
