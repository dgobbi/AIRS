/*=========================================================================

  Program:   Atamai Image Registration and Segmentation
  Module:    lsmtonifti.cxx

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


#include "vtkNIFTIWriter.h"
#include "vtkLSMReader.h"

#include <vtkMath.h>
#include <vtkImageClip.h>
#include <vtkImageChangeInformation.h>
#include <vtkImageResize.h>
#include <vtkImageHistogramStatistics.h>
#include <vtkImageGridSource.h>
#include <vtkImageAppend.h>
#include <vtkImageAppendComponents.h>
#include <vtkImageBlend.h>
#include <vtkImageShiftScale.h>
#include <vtkImageMask.h>
#include <vtkImageThreshold.h>
#include <vtkMatrix4x4.h>
#include <vtkImageData.h>
#include <vtkIdTypeArray.h>
#include <vtkErrorCode.h>
#include <vtkSmartPointer.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <ctype.h>

#define LSMTONIFTI_VERSION "1.0.0"

// Get the command name without the full path
const char *basename(const char *command_name)
{
  const char *cp = command_name + strlen(command_name);
  while (cp != command_name && cp[-1] != '\\' && cp[-1] != '/')
    {
    --cp;
    }
  return cp;
}

// Print the options
void usage(FILE *file, const char *command_name)
{
  fprintf(file,
    "%s version %s\n"
    "Convert Zeiss .lsm files into NIfTI .nii files.\n"
    "usage: %s [options] input2.lsm [output.nii[.gz]]\n",
    command_name, LSMTONIFTI_VERSION, command_name);
}

// Print the version
void version(FILE *file)
{
  fprintf(file, "%s\n", LSMTONIFTI_VERSION);
}

// Print a short help message
void shorthelp(FILE *file, const char *command_name)
{
  usage(file, command_name);

  fprintf(file, "\n");
  fprintf(file,
    "Options:\n"
    "  --normalize      Normalize the slice intensity while converting.\n"
    "  --crop WxH+X+Y   Specify a cropping rectangle.\n"
    "  --resize WxH     Resize the image to the specified dimensions.\n"
    "  --grid WxH[+X+Y] Overlay a grid with spacing WxH and offset X,Y.\n"
    "  --version        Print the version and exit.\n"
    "  --help           Print documentation.\n\n"
  );
}

// Print additional help
void longhelp(FILE *file, const char *command_name)
{
  usage(file, command_name);

  fprintf(file, "\n");
  fprintf(file,
    "This program converts .lsm files into NIfTI format.  Please note that only\n"
    "uncompressed .lsm files are supported.  If no output filename is specified,\n"
    "then one will be created with the same name as the input file, but with the\n"
    "original extension replaced by a \".nii\" extension.\n\n"
    "Several options are available for the conversion:\n\n");
  fprintf(file,
    "--normalize        Normalize the intensity of each slice to match that of\n"
    "                   the middle slice in the stack.  The intensities will be\n"
    "                   scaled to optimize the similarity, in a least-squares\n"
    "                   sense, of the slices to each other.\n\n");
  fprintf(file,
    "--crop WxH+X+Y     Crop the image to a WxH region with an upper left corner\n"
    "                   at X,Y.  If X=0 and Y=0, or if +X+Y is omitted, then the\n"
    "                   upper left corner of the original image will be used.\n\n");
  fprintf(file,
    "--resize WxH       Resize the image to WxH.  The image will be bandpass-\n"
    "                   filtered to avoid aliasing and resampled with a Lanczos\n"
    "                   interpolation kernel.\n\n");
  fprintf(file,
    "--grid WxH[+X+Y]   Overlay a light-grey grid on the image where the size of\n"
    "                   each square in the grid is WxH.  Optionally, an offset of\n"
    "                   X,Y can also be specified.\n\n");
}

// Parse a geometry option in the format "WxH+X+Y"
void geometry(char *argv[], int argc, int argi, int g[4])
{
  const char *option = argv[argi++];

  if (argi == argc || argv[argi][0] == '-')
    {
    const char *command_name = basename(argv[0]);
    fprintf(stderr, "%s: option %s must be followed by a value.\n",
            command_name, option);
    shorthelp(stderr, command_name);
    exit(1);
    }

  char *arg = argv[argi];
  const char delim[4] = { 'x', '+', '+', '\0' };
  for (int i = 0; i < 4; i++)
    {
    g[i] = static_cast<int>(strtoul(arg, &arg, 10));
    if (*arg == delim[i])
      {
      arg++;
      }
    else if (i == 1 && *arg == '\0')
      {
      g[2] = 0;
      g[3] = 0;
      break;
      }
    else
      {
      const char *command_name = basename(argv[0]);
      fprintf(stderr, "%s: option %s must be followed by a valid geometry.\n",
              command_name, option);
      shorthelp(stderr, command_name);
      exit(1);
      }
    }
}

// Check for IO errors that might have occurred in VTK classes
void check_error(vtkObject *o)
{
  vtkLSMReader *reader = vtkLSMReader::SafeDownCast(o);
  vtkNIFTIWriter *writer = vtkNIFTIWriter::SafeDownCast(o);
  const char *filename = 0;
  unsigned long errorcode = 0;

  if (writer)
    {
    filename = writer->GetFileName();
    errorcode = writer->GetErrorCode();
    }
  else if (reader)
    {
    filename = reader->GetFileName();
    errorcode = reader->GetErrorCode();
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

// Compute the squared difference between two histograms,
// ignore first and last bins because these have saturation
double histdiff(const int hist1[], const int hist2[], int n, double scale)
{
  double s = 0.0;
  double c = 0.0;

  for (int i = 1; i < n-1; i++)
    {
    int j = static_cast<int>(i/scale);
    if (j > 0 && j < n-1)
      {
      double d = hist1[i]*scale - hist2[j];
      s += d*d;
      c += 1.0;
      }
    }

  return s/c;
}

// Match two histograms, return the scale factor
double histmatch(const int hist1[], const int hist2[], int n)
{
  // The search is done on a logarithmic scale, rather than
  // a linear scale.

  // search range for scale is [0.1, 10.0]
  double erange = log(10.0);
  // initial guess for scale is 1.0
  double ebase = 0.0;
  double bestscale = 1.0;
  double bestsum = 1e30;
  // at each level, try 101 values evenly spaced on the log scale
  int steps = 101;
  double halfwidth = 0.5*(steps - 1);

  // start by searching the full range, and reduce the range each iteration
  for (int i = 0; i < 10; i++)
    {
    // find the best value within the current range
    for (int ei = 0; ei < steps; ei++)
      {
      // e will go from ebase-erange to ebase+erange
      double e = ebase + (ei - halfwidth)/halfwidth*erange;
      double scale = exp(e);
      double s = histdiff(hist1, hist2, n, scale);
      if (s < bestsum)
        {
        bestsum = s;
        bestscale = scale;
        }
      }
    // use ebase as centre of search range for the next time around
    ebase = log(bestscale);
    // reduce the width of the search range for the next time around
    erange /= halfwidth;
    }

  return bestscale;
}

// Main program
int main(int argc, char *argv[])
{
  // get the command name
  const char *command_name = basename(argv[0]);

  // check options
  bool normalizeOption = false;
  bool cropOption = false;
  int cropGeometry[4] = { 1024, 1024, 0, 0 };
  bool resizeOption = false;
  int resizeGeometry[4] = { 1024, 1024, 0, 0 };
  bool gridOption = false;
  int gridGeometry[4] = { 10, 10, 0, 0 };
  bool compressOption = false;

  const char *infile = 0;
  const char *outfile = 0;

  int argi = 1;
  bool optflag = true;
  for (; argi < argc; argi++)
    {
    const char *arg = argv[argi];

    if (optflag && arg[0] == '-')
      {
      if (strcmp(arg, "--") == 0)
        {
        // treat all remaining options as filenames
        optflag = false;
        }
      else if (strcmp(arg, "-z") == 0)
        {
        compressOption = true;
        }
      else if (strcmp(arg, "--normalize") == 0)
        {
        normalizeOption = true;
        }
      else if (strcmp(arg, "--crop") == 0)
        {
        cropOption = true;
        geometry(argv, argc, argi++, cropGeometry);
        }
      else if (strcmp(arg, "--resize") == 0)
        {
        resizeOption = true;
        geometry(argv, argc, argi++, resizeGeometry);
        }
      else if (strcmp(arg, "--grid") == 0)
        {
        gridOption = true;
        geometry(argv, argc, argi++, gridGeometry);
        }
      else if (strcmp(arg, "--help") == 0)
        {
        longhelp(stdout, command_name);
        exit(0);
        }
      else if (strcmp(arg, "-h") == 0)
        {
        shorthelp(stdout, command_name);
        exit(0);
        }
      else if (strcmp(arg, "--version") == 0)
        {
        version(stdout);
        exit(0);
        }
      else
        {
        fprintf(stderr, "%s: unrecognized option %s.\n", command_name, arg);
        shorthelp(stderr, command_name);
        exit(1);
        }
      }
    else
      {
      if (infile == 0)
        {
        infile = arg;
        }
      else if (outfile == 0)
        {
        outfile = arg;
        }
      else
        {
        fprintf(stderr, "%s: Too many arguments.\n", command_name);
        usage(stderr, command_name);
        exit(1);
        }
      }
    }

  // check for an input file name
  if (infile == 0)
    {
    fprintf(stderr, "%s: No input file was provided.\n", command_name);
    usage(stderr, command_name);
    exit(1);
    }

  // generate an output file name if not provided
  char *outfile_generated = 0;
  if (outfile == 0)
    {
    const char *ext2[2] = { ".nii", ".nii.gz" };
    const char *ext = ext2[compressOption];

    size_t l = strlen(infile);
    size_t n = l;
    while (n > 0)
      {
      --n;
      if (infile[n] == '.')
        {
        break;
        }
      else if (!isalnum(infile[n]))
        {
        n = 0;
        }
      }
    if (n == 0)
      {
      n = l;
      }

    size_t m = strlen(ext);
    outfile_generated = new char [n + m + 1];
    strncpy(outfile_generated, infile, n);
    strcpy(&outfile_generated[n], ext);
    outfile = outfile_generated;
    }
  else
    {
    size_t l = strlen(outfile);
    if (l > 3 && strcmp(outfile+(l-3), ".gz") == 0)
      {
      l -= 3;
      }
    if (l < 4 || strncmp(outfile+(l-4), ".nii", 4) != 0)
      {
      fprintf(stderr, "%s: Output filename must end in .nii or .nii.gz.\n",
              command_name);
      usage(stderr, command_name);
      exit(1);
      }
    }

  // for connecting to next filter
  vtkAlgorithmOutput *port = 0;

  // read the lsm file
  vtkSmartPointer<vtkLSMReader> reader =
    vtkSmartPointer<vtkLSMReader>::New();
  reader->SetFileName(infile);
  reader->Update();
  check_error(reader);
  port = reader->GetOutputPort();

  // get the dimensions of the image
  int extent[6];
  reader->GetOutput()->GetExtent(extent);

  // check if the user requested for the image to be cropped
  vtkSmartPointer<vtkImageClip> cropper =
    vtkSmartPointer<vtkImageClip>::New();
  if (cropOption)
    {
    if (cropGeometry[0] == 0 ||
        cropGeometry[1] == 0 ||
        cropGeometry[0] + cropGeometry[2] > extent[1] - extent[0] + 1 ||
        cropGeometry[1] + cropGeometry[3] > extent[3] - extent[2] + 1)
      {
      fprintf(stderr, "%s: Crop region is outside of image bounds for %s.\n",
              command_name, infile);
      exit(1);
      }

    extent[0] += cropGeometry[2];
    extent[1] = extent[0] + cropGeometry[0] + 1;
    extent[2] += cropGeometry[3];
    extent[3] = extent[2] + cropGeometry[1] + 1;

    // crop the image
    cropper->SetInputConnection(port);
    cropper->SetOutputWholeExtent(extent);
    cropper->ClipDataOn();
    cropper->Update();
    port = cropper->GetOutputPort();
    }

  // shift the extent so that it starts at (0,0,0)
  vtkSmartPointer<vtkImageChangeInformation> changer =
    vtkSmartPointer<vtkImageChangeInformation>::New();
  changer->SetInputConnection(port);
  changer->SetOutputOrigin(0.0, 0.0, 0.0);
  changer->SetOutputExtentStart(0, 0, 0);
  changer->Update();
  port = changer->GetOutputPort();

  // get the new, shifted extent
  changer->GetOutput()->GetExtent(extent);

  // perform the intensity normalization
  vtkSmartPointer<vtkImageAppend> append =
    vtkSmartPointer<vtkImageAppend>::New();
  append->SetAppendAxis(2);
  if (normalizeOption)
    {
    // compute an intensity scaling factor for each slice
    int numslices = extent[5] - extent[4] + 1;
    double *scales = new double[numslices];

    // build a histogram for each slice
    int **allhists = new int *[numslices];
    int histsize = 0;

    for (int i = extent[4]; i <= extent[5]; i++)
      {
      vtkSmartPointer<vtkImageClip> clip =
        vtkSmartPointer<vtkImageClip>::New();
      clip->SetInputConnection(port);
      clip->SetOutputWholeExtent(extent[0], extent[1],
                                 extent[2], extent[3],
                                 i, i);
      clip->ClipDataOn();

      vtkSmartPointer<vtkImageHistogramStatistics> hist =
        vtkSmartPointer<vtkImageHistogramStatistics>::New();
      hist->SetInputConnection(clip->GetOutputPort());
      hist->Update();
      vtkIdTypeArray *histdata = hist->GetHistogram();

      histsize = static_cast<int>(histdata->GetNumberOfTuples());
      int *h = new int[histsize];
      for (int j = 0; j < histsize; j++)
        {
        h[j] = static_cast<int>(histdata->GetValue(j));
        }

      allhists[i - extent[4]] = h;
      }

    // compute the scaling for each histogram
    // (match them all to the central slice)
    int midhist = numslices/2;
    for (int k = 0; k < numslices; k++)
      {
      scales[k] = histmatch(allhists[k], allhists[midhist], histsize);
      }

    for (int k = 0; k < numslices; k++)
      {
      delete [] allhists[k];
      }
    delete [] allhists;

    // correct the intensity range of each slice
    for (int i = extent[4]; i <= extent[5]; i++)
      {
      vtkSmartPointer<vtkImageClip> clip =
        vtkSmartPointer<vtkImageClip>::New();
      clip->SetInputConnection(port);
      clip->SetOutputWholeExtent(extent[0], extent[1],
                                 extent[2], extent[3],
                                 i, i);
      clip->ClipDataOn();

      vtkSmartPointer<vtkImageShiftScale> rescale =
        vtkSmartPointer<vtkImageShiftScale>::New();
      rescale->SetInputConnection(clip->GetOutputPort());
      rescale->SetScale(1.0/scales[i - extent[4]]);
      rescale->ClampOverflowOn();
      rescale->Update();

      // leave saturated pixels at the highest possible value
      // when rescaling the intensity?  off by default.
      bool maintain_saturation = false;
      if (maintain_saturation)
        {
        vtkSmartPointer<vtkImageThreshold> thresh =
          vtkSmartPointer<vtkImageThreshold>::New();
        thresh->SetInputConnection(clip->GetOutputPort());
        thresh->ThresholdByUpper(255);
        thresh->ReplaceOutOn();
        thresh->SetOutValue(0);
        thresh->Update();

        vtkSmartPointer<vtkImageMask> mask =
          vtkSmartPointer<vtkImageMask>::New();
        mask->SetInputConnection(0, rescale->GetOutputPort());
        mask->SetInputConnection(1, thresh->GetOutputPort());
        mask->SetMaskedOutputValue(255);
        mask->NotMaskOn();
        mask->Update();

        append->AddInputConnection(mask->GetOutputPort());
        }
      else
        {
        append->AddInputConnection(rescale->GetOutputPort());
        }
      }

    port = append->GetOutputPort();

    delete [] scales;
    }

  // resize the image
  vtkSmartPointer<vtkImageResize> resize =
    vtkSmartPointer<vtkImageResize>::New();
  if (resizeOption)
    {
    resize->SetInputConnection(port);
    resize->SetBorder(1);
    resize->SetResizeMethodToOutputDimensions();
    resize->SetOutputDimensions(resizeGeometry[0], resizeGeometry[1],
                                (extent[5] - extent[4] + 1));
    resize->Update();
    port = resize->GetOutputPort();
    }

  // create the grid overlay
  vtkSmartPointer<vtkImageBlend> blend =
    vtkSmartPointer<vtkImageBlend>::New();
  if (gridOption)
    {
    vtkSmartPointer<vtkImageChangeInformation> info =
      vtkSmartPointer<vtkImageChangeInformation>::New();
    info->SetInputConnection(port);
    info->Update();

    vtkSmartPointer<vtkImageGridSource> grid =
      vtkSmartPointer<vtkImageGridSource>::New();
    grid->SetGridSpacing(gridGeometry[0], gridGeometry[1], 0);
    grid->SetGridOrigin(gridGeometry[2], gridGeometry[3], 0);
    grid->SetDataOrigin(info->GetOutput()->GetOrigin());
    grid->SetDataSpacing(info->GetOutput()->GetSpacing());
    grid->SetDataExtent(info->GetOutput()->GetExtent());
    grid->SetLineValue(255);
    grid->SetFillValue(0);
    grid->SetDataScalarTypeToUnsignedChar();

    vtkSmartPointer<vtkImageAppendComponents> appendc =
      vtkSmartPointer<vtkImageAppendComponents>::New();
    appendc->SetInputConnection(grid->GetOutputPort());
    appendc->AddInputConnection(grid->GetOutputPort());

    blend->SetInputConnection(port);
    blend->AddInputConnection(appendc->GetOutputPort());
    blend->SetOpacity(1, 0.5);
    blend->Update();
    port = blend->GetOutputPort();
    }

  // create an identity matrix for the NIfTI header
  vtkSmartPointer<vtkMatrix4x4> matrix =
    vtkSmartPointer<vtkMatrix4x4>::New();

  // write the nifti file
  vtkSmartPointer<vtkNIFTIWriter> writer =
    vtkSmartPointer<vtkNIFTIWriter>::New();
  writer->SetInputConnection(port);
  writer->SetFileName(outfile);
  writer->SetSFormMatrix(matrix);
  writer->SetQFormMatrix(matrix);
  writer->Write();
  check_error(reader);

  delete [] outfile_generated;

  return 0;
}
