/*=========================================================================

  Program:   Atamai Image Registration and Segmentation
  Module:    WarpImage.cxx

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

#include <vtkImageData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkImageReslice.h>
#include <vtkMNITransformReader.h>
#include <vtkGeneralTransform.h>
#include <vtkGridTransform.h>
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
void usage(FILE *file, const char *command_name)
{
  const char *cp = command_name + strlen(command_name);
  while (cp != command_name && cp[-1] != '\\' && cp[-1] != '/') { --cp; }

  fprintf(file,
    "usage: %s Input.nii -o Output.nii -R Like.nii Warp.nii Affine.txt\n", cp);
  fprintf(file,
    "usage: %s Input.nii -o Output.nii -R Like.nii -i Affine.txt InverseWarp.nii\n", cp);
  fprintf(file,
    "options:\n"
    "  --version               Print the version and exit.\n"
    "  --help                  Print minimal documentation.\n"
  );
}

// Print the help
void help(FILE *file, const char *command_name)
{
  usage(file, command_name);

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
void check_error(vtkObject *o)
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

int read_transform(
  vtkGeneralTransform *transform, const char *file, bool invert)
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
    check_error(reader);
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
    check_error(reader);
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
    check_error(reader);

    // break the pipeline connection to the reader
    vtkSmartPointer<vtkImageData> image =
      vtkSmartPointer<vtkImageData>::New();
    image->CopyStructure(reader->GetOutput());
    image->GetPointData()->PassData(reader->GetOutput()->GetPointData());

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

    vtkSmartPointer<vtkGridTransform> gt =
      vtkSmartPointer<vtkGridTransform>::New();
    // use linear to match ANTS?
    gt->SetInterpolationModeToLinear();
    gt->SetDisplacementGrid(image);
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

  // the output and target file names
  const char *infile = 0;
  const char *outfile = 0;
  const char *targetfile = 0;

  // read the options from the command line
  int argi = 1;
  while (argi < argc)
    {
    const char *arg = argv[argi++];
    if (strcmp(arg, "--help") == 0)
      {
      help(stdout, argv[0]);
      exit(0);
      }
    else if (strcmp(arg, "-o") == 0)
      {
      if (argi >= argc)
        {
        fprintf(stderr, "\nA file must follow the \'-o\' flag\n\n");
        usage(stderr, argv[0]);
        exit(1);
        }
      arg = argv[argi++];
      outfile = arg;
      }
    else if (strcmp(arg, "-R") == 0)
      {
      if (argi >= argc)
        {
        fprintf(stderr, "\nA file must follow the \'-R\' flag\n\n");
        usage(stderr, argv[0]);
        exit(1);
        }
      arg = argv[argi++];
      targetfile = arg;
      }
    else if (strcmp(arg, "-i") == 0)
      {
      if (argi >= argc)
        {
        fprintf(stderr, "\nA file must follow the \'-R\' flag\n\n");
        usage(stderr, argv[0]);
        exit(1);
        }
      arg = argv[argi++];
      if (!read_transform(transform, arg, true))
        {
        exit(1);
        }
      }
    else if (arg[0] != '-')
      {
      if (!infile)
        {
        infile = arg;
        }
      else if (!read_transform(transform, arg, false))
        {
        exit(1);
        }
      }
    else
      {
      fprintf(stderr, "\nUnrecognized option %s.\n\n", arg);
      usage(stderr, argv[0]);
      exit(1);
      }
    }

  if (!infile)
    {
    fprintf(stderr, "\nAn input file must be specified.\n");
    exit(1);
    }

  if (!outfile)
    {
    fprintf(stderr, "\nAn output file must be specified with \"-o\"\n");
    exit(1);
    }

  vtkSmartPointer<vtkNIFTIReader> reader =
    vtkSmartPointer<vtkNIFTIReader>::New();
  reader->SetFileName(infile);
  reader->Update();
  check_error(reader);

  vtkSmartPointer<vtkNIFTIReader> rreader =
    vtkSmartPointer<vtkNIFTIReader>::New();
  if (targetfile)
    {
    rreader->SetFileName(targetfile);
    rreader->Update();
    check_error(rreader);
    }
  else
    {
    rreader = reader;
    }

  double spacing[3], origin[3];
  int extent[6];
  vtkImageData *image = rreader->GetOutput();
  image->GetSpacing(spacing);
  image->GetOrigin(origin);
  image->GetExtent(extent);

  vtkSmartPointer<vtkImageReslice> reslice =
    vtkSmartPointer<vtkImageReslice>::New();
  reslice->SetInputConnection(reader->GetOutputPort());
  reslice->SetResliceTransform(transform);
  reslice->SetOutputSpacing(spacing);
  reslice->SetOutputOrigin(origin);
  reslice->SetOutputExtent(extent);
  reslice->SetInterpolationModeToLinear();
  reslice->Update();

  vtkSmartPointer<vtkNIFTIWriter> writer =
    vtkSmartPointer<vtkNIFTIWriter>::New();
  writer->SetInputConnection(reslice->GetOutputPort());
  writer->SetFileName(outfile);
  writer->Write();
  check_error(writer);
}
