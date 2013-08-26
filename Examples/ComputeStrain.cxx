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

#include <vtkMatrix4x4.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkStringArray.h>
#include <vtkIntArray.h>
#include <vtkInformation.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTransform.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageExtractComponents.h>
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
    "usage: %s -o Output.nii [options] Warp.nii Affine.txt [...]\n", cp);
  fprintf(file,
    "usage: %s -o Output.nii [options] -i Affine.txt InverseWarp.nii [...]\n", cp);
  fprintf(file,
    "\noptions:\n"
    "  -o <output.nii[.gz]>    The name of the output file.\n"
    "  -R <file.nii[.gz]>      A file with the desired output size and spacing.\n"
    "  -i Affine.txt           Invert the affine transformation that follows.\n"
    "  --size WxH[xD]          Specify the desired size of the output file.\n"
    "  --deformation-gradient  Output the deformation gradient tensor.\n"
    "  --greens-strain-tensor  Output Green's strain tensor (the default).\n"
    "  --principal-directions  Output the principal directions of Green's strain.\n"
    "  --principal-components  Output the principal components of Green's straing.\n"
    "  --principal-component   Output just the largest principal component.\n"
    "  --help                  Print a short help document.\n"
  );
}

// Print the help
void strain_help(FILE *file, const char *command_name)
{
  strain_usage(file, command_name);

  fprintf(file,
    "\n");

  fprintf(file,
    "This program computes strain tensors from a series of deformations.\n"
    "\n");
  fprintf(file,
    "The deformations are described in two ways: as an affine transformation\n"
    "stored as a 3x3 matrix and an offset vector, or as a displacement vector\n"
    "field stored in a nifti (.nii) file.  When the displacement vector field is\n"
    "loaded into this program from the nifti file, a three-dimensional cubic\n"
    "B-spline is computed for the x, y, and z components of the displacement\n"
    "in order to create a smooth, differentiable vector field.\n"
    "\n");
  fprintf(file,
    "In general, the full deformation of an image is composed of a displacement\n"
    "vector field that describes the fine structure of the deformation,\n"
    "followed by an affine transformation that provides the global initial\n"
    "estimate of the deformation.  This order is used, rather than the reverse,\n"
    "because it allows the samples in the vector field to be aligned with the\n"
    "image that is being deformed.\n"
    "\n");
  fprintf(file,
    "It is also possible to describe a deformation as a concatenation of\n"
    "multiple affine transformations and deformation vector fields.  To do\n"
    "this, simply supply all the transformation files to be concatenated on\n"
    "the command line.  The strain tensor will be computed from the total\n"
    "deformation that results from applying each sub-deformation in turn.\n"
    "\n");
  fprintf(file,
    "The --size option sets the size of the output image to WxHxD.  If the\n"
    "depth (xD) is omitted, then the number of slices in the output tensor image\n"
    "will be the same as the number of slices in the original image stack.  If\n"
    "the requested resolution is coarser (i.e. the WxHxD size is smaller) than\n"
    "a loaded displacement field files, then a guassian smoothing will be applied\n"
    "to the displacement field before it is used to compute the tensors.  The\n"
    "sigma_x, sigma_y, and sigma_z of the Gaussian used for the smoothing will\n"
    "equal to 0.399 times the ratio of the output sample spacing (tensor spacing)\n"
    "to the input sample spacing (displacement vector spacing).  This smoothing\n"
    "eliminates the aliasing that would otherwise result from undersampling the\n"
    "displacement vector field while computing the tensors.\n"
    "\n");
  fprintf(file,
    "If an image file is supplied with the -R option, then the strain will be\n"
    "computed at every pixel in this supplied image.  For example, if you supply\n"
    "one of the original microscope image stacks with -R, then the output tensor\n"
    "image will have the same geometry as that file.\n"
    "\n");
  fprintf(file,
    "Four output options are available:\n\n");
  fprintf(file,
    "--deformation-gradient computes a 3x3 matrix containing the partial\n"
    "derivatives of the deformed coordinates (x\',y\',z\') with respect to\n"
    "the original, undeformed coordinates (x,y,z).\n\n");
  fprintf(file,
    "--greens-strain-tensor computes Green\'s strain tensor via the\n"
    "deformation gradient: E = 0.5*(F\'F - I) where F is the deformation\n"
    "gradient tensor and I is the 3x3 identity matrix.\n\n");
  fprintf(file,
    "--principal-directions computes the principal directions of Green\'s\n"
    "strain tensor via the Jacobi algorithm.\n\n");
  fprintf(file,
    "--principal-components computes the principal components of Green\'s\n"
    "strain tensor, and outputs them in order from largest to smallest.\n"
    "The order will match the principal directions/\n\n");
  fprintf(file,
    "--principal-component will output only the largest principal component.\n"
    "\n");
  fprintf(file,
    "The name of the desired output file must be preceded by \"-o\".\n"
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

// Parse a geometry option in the format "WxHxD"
void dimensions(char *argv[], int argc, int argi, int g[3])
{
  const char *option = argv[argi++];

  if (argi == argc || argv[argi][0] == '-')
    {
    fprintf(stderr, "\nError: option %s must be followed by a value.\n",
            option);
    strain_usage(stderr, argv[0]);
    exit(1);
    }

  char *arg = argv[argi];
  for (int i = 0; i < 3; i++)
    {
    g[i] = static_cast<int>(strtoul(arg, &arg, 10));
    if (*arg == 'x')
      {
      arg++;
      }
    else if (i == 1 && *arg == '\0')
      {
      g[2] = -1;
      break;
      }
    else
      {
      fprintf(stderr,
              "\nError: option %s must be followed by valid dimensions.\n",
              option);
      strain_usage(stderr, argv[0]);
      exit(1);
      }
    }
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
  int outputSize[3] = { -1, -1, -1 };

  enum OutputType {
    DeformationGradient,
    GreensStrainTensor,
    PrincipalDirections,
    PrincipalComponents,
    PrincipalComponent
  };
  OutputType outputType = GreensStrainTensor;

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
    else if (strcmp(arg, "--deformation-gradient") == 0)
      {
      outputType = DeformationGradient;
      }
    else if (strcmp(arg, "--greens-strain-tensor") == 0)
      {
      outputType = GreensStrainTensor;
      }
    else if (strcmp(arg, "--principal-directions") == 0)
      {
      outputType = PrincipalDirections;
      }
    else if (strcmp(arg, "--principal-components") == 0)
      {
      outputType = PrincipalComponents;
      }
    else if (strcmp(arg, "--principal-component") == 0)
      {
      outputType = PrincipalComponent;
      }
    else if (strcmp(arg, "--size") == 0)
      {
      dimensions(argv, argc, argi++, outputSize);
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
    // if no targetfile, use the first deformation field file
    if (transformFiles->GetNumberOfValues() > 0)
      {
      const char *firstTrans = transformFiles->GetValue(0);
      size_t l = strlen(firstTrans);
      if ((l >= 4 && strcmp(&firstTrans[l-4], ".nii") == 0) ||
          (l >= 7 && strcmp(&firstTrans[l-7], ".nii.gz") == 0))
        {
        targetfile = firstTrans;
        }
      }
    if (!targetfile)
      {
      fprintf(stderr, "\nA target file must be specified with \"-R\" "
              "unless the first transform file is a Warp.nii file.\n");
      exit(1);
      }
    }

  vtkSmartPointer<vtkNIFTIReader> reader =
    vtkSmartPointer<vtkNIFTIReader>::New();
  reader->SetFileName(targetfile);
  reader->UpdateInformation();
  strain_check_error(reader);

  double spacing[3], origin[3];
  int extent[6];
  vtkInformation *info = reader->GetOutputPortInformation(0);
  info->Get(vtkDataObject::SPACING(), spacing);
  info->Get(vtkDataObject::ORIGIN(), origin);
  info->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent);

  // adjust the sample spacing to get the desired output image size
  for (int j = 0; j < 3; j++)
    {
    if (outputSize[j] >= 0)
      {
      double o = origin[j] + spacing[j]*(extent[2*j] - 0.5);
      double d = (extent[2*j + 1] - extent[2*j] + 1)*spacing[j];
      spacing[j] = d/outputSize[j];
      extent[2*j + 1] = extent[2*j] + outputSize[j] - 1;
      origin[j] = o - spacing[j]*(extent[2*j] - 0.5);
      }
    }

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
  switch (outputType)
    {
    case DeformationGradient:
      computeStrain->SetOutputValueToDeformationGradient();
      break;
    case GreensStrainTensor:
      computeStrain->SetOutputValueToGreensStrainTensor();
      break;
    case PrincipalDirections:
      computeStrain->SetOutputValueToPrincipalDirections();
      break;
    case PrincipalComponents:
    case PrincipalComponent:
      computeStrain->SetOutputValueToPrincipalComponents();
      break;
    }
  computeStrain->SetOutputScalarTypeToFloat();
  computeStrain->SetInput(transform);
  computeStrain->SetOutputSpacing(spacing);
  computeStrain->SetOutputOrigin(origin);
  computeStrain->SetOutputExtent(extent);
  computeStrain->Update();

  vtkSmartPointer<vtkImageExtractComponents> extractor =
    vtkSmartPointer<vtkImageExtractComponents>::New();
  extractor->SetInput(computeStrain->GetOutput());
  if (outputType == PrincipalComponents)
    {
    extractor->SetComponents(0,1,2);
    }
  else if (outputType == PrincipalComponent)
    {
    extractor->SetComponents(0);
    }

  vtkSmartPointer<vtkNIFTIWriter> writer =
    vtkSmartPointer<vtkNIFTIWriter>::New();
  if (outputType == PrincipalComponents ||
      outputType == PrincipalComponent)
    {
    writer->SetInputConnection(extractor->GetOutputPort());
    }
  else
    {
    writer->SetInputConnection(computeStrain->GetOutputPort());
    }
  writer->SetFileName(outfile);
  writer->Write();
  strain_check_error(writer);
}
