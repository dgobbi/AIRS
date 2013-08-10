/*=========================================================================

Program:   Atamai Image Registration and Segmentation
Module:    register.cxx

   This software is distributed WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/

// This example demonstrates rigid registration of images.  Since it is a
// rigid registration, the expectation is that the images come from the
// same patient.

// Two image file formats are supported for this example: MINC and DICOM.
// DICOM images are read with the troublesome vtkDICOMImageReader, which
// may get the slice spacing or ordering wrong, or even fail to read the
// images altogether.

// Image registration is done first on a blurred, low-resolution version of
// the image before being done on the full resolution image, and is also
// done first with no interpolation before being done with linear interpolation.
// This multi-stage approach increases the robustness and often the speed of
// the registration.

#include <vtkSmartPointer.h>

#include <vtkImageReslice.h>
#include <vtkImageResize.h>
#include <vtkImageSincInterpolator.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkMath.h>

#include <vtkMINCImageReader.h>
#include <vtkMNITransformReader.h>
#include <vtkMNITransformWriter.h>
#include <vtkDICOMImageReader.h>

#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageSlice.h>
#include <vtkImageStack.h>
#include <vtkImageResliceMapper.h>
#include <vtkImageProperty.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkJPEGWriter.h>

#include <vtkTimerLog.h>

#include "AIRSConfig.h"
#include "vtkITKXFMReader.h"
#include "vtkITKXFMWriter.h"
#include "vtkImageRegistration.h"

// optional readers
#ifdef AIRS_USE_DICOM
#define AIRS_USE_NIFTI
#include <vtkNIFTIReader.h>
#include <vtkDICOMReader.h>
#include <vtkDICOMSorter.h>
#include <vtkGlobFileNames.h>
#endif

// coord systems
enum { NativeCoords, DICOMCoords, NIFTICoords };

// internal methods for reading images, these methods read the image
// into the specified data object and also provide a matrix for converting
// the data coordinates into patient coordinates.
namespace {

#ifdef AIRS_USE_DICOM
void ReadDICOMImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *directoryName,
  int coordSystem)
{
  // get the files
  vtkSmartPointer<vtkGlobFileNames> glob =
    vtkSmartPointer<vtkGlobFileNames>::New();
  glob->SetDirectory(directoryName);
  glob->AddFileNames("*");

  // sort the files
  vtkSmartPointer<vtkDICOMSorter> sorter =
    vtkSmartPointer<vtkDICOMSorter>::New();
  sorter->SetInputFileNames(glob->GetFileNames());
  sorter->Update();

  if (sorter->GetNumberOfSeries() != 1)
    {
    fprintf(stderr, "directory contains %d DICOM series, expected 1: %s\n",
            static_cast<int>(sorter->GetNumberOfSeries()), directoryName);
    exit(1);
    }

  // read the image
  vtkSmartPointer<vtkDICOMReader> reader =
    vtkSmartPointer<vtkDICOMReader>::New();
  reader->SetFileNames(sorter->GetFileNamesForSeries(0));
  if (coordSystem == NIFTICoords)
    {
    reader->SetMemoryRowOrderToBottomUp();
    }
  else
    {
    reader->SetMemoryRowOrderToFileNative();
    }
  reader->Update();

  vtkImageData *image = reader->GetOutput();

  // get the data
  data->CopyStructure(image);
  data->GetPointData()->PassData(image->GetPointData());

  // get the matrix
  matrix->DeepCopy(reader->GetPatientMatrix());
}

#else

void ReadDICOMImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *directoryName,
  int coordSystem)
{
  // read the image
  vtkSmartPointer<vtkDICOMImageReader> reader =
    vtkSmartPointer<vtkDICOMImageReader>::New();

  reader->SetDirectoryName(directoryName);
  reader->Update();

  vtkSmartPointer<vtkImageData> image = reader->GetOutput();

  if (coordSystem != NIFTICoords)
    {
    // the reader flips the image and reverses the ordering, so undo these
    vtkSmartPointer<vtkImageReslice> flip =
      vtkSmartPointer<vtkImageReslice>::New();

    flip->SetInputConnection(reader->GetOutputPort());
    flip->SetResliceAxesDirectionCosines(
      1,0,0, 0,-1,0, 0,0,-1);
    flip->Update();

    image = flip->GetOutput();
    }

  // get the data
  data->CopyStructure(image);
  data->GetPointData()->PassData(image->GetPointData());
  data->SetOrigin(0,0,0);

  // generate the matrix
  float *position = reader->GetImagePositionPatient();
  float *orientation = reader->GetImageOrientationPatient();
  float *xdir = &orientation[0];
  float *ydir = &orientation[3];
  float zdir[3];
  vtkMath::Cross(xdir, ydir, zdir);

  for (int i = 0; i < 3; i++)
    {
    matrix->Element[i][0] = xdir[i];
    matrix->Element[i][1] = ydir[i];
    matrix->Element[i][2] = zdir[i];
    matrix->Element[i][3] = position[i];
    }
  matrix->Element[3][0] = 0;
  matrix->Element[3][1] = 0;
  matrix->Element[3][2] = 0;
  matrix->Element[3][3] = 1;

  if (coordSystem == NIFTICoords)
    {
    double spacing[3], origin[3];
    int extent[6];
    image->GetSpacing(spacing);
    image->GetOrigin(origin);
    image->GetExtent(extent);
    // account fo the y and z flips
    double point[4];
    point[0] = origin[0] + spacing[0]*extent[0];
    point[1] = origin[1] + spacing[1]*extent[3];
    point[2] = origin[2] + spacing[2]*extent[5];
    point[3] = 1.0;
    matrix->MultiplyPoint(point, point);
    for (int j = 0; j < 3; j++)
      {
      matrix->Element[j][1] = -matrix->Element[j][1];
      matrix->Element[j][2] = -matrix->Element[j][2];
      matrix->Element[j][3] = point[j];
      }
    // do the DICOM to NIFTI coord conversion
    for (int k = 0; k < 4; k++)
      {
      matrix->Element[0][k] = -matrix->Element[0][k];
      matrix->Element[1][k] = -matrix->Element[1][k];
      }
    }

  matrix->Modified();
}
#endif

void ReadMINCImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *fileName,
  int coordSystem)
{
  // read the image
  vtkSmartPointer<vtkMINCImageReader> reader =
    vtkSmartPointer<vtkMINCImageReader>::New();

  reader->SetFileName(fileName);
  reader->Update();

  vtkSmartPointer<vtkImageData> image = reader->GetOutput();

  if (coordSystem == DICOMCoords)
    {
    double spacing[3];
    reader->GetOutput()->GetSpacing(spacing);
    spacing[0] = fabs(spacing[0]);
    spacing[1] = fabs(spacing[1]);
    spacing[2] = fabs(spacing[2]);

    // flip the image rows into a DICOM-style ordering
    vtkSmartPointer<vtkImageReslice> flip =
      vtkSmartPointer<vtkImageReslice>::New();

    flip->SetInputConnection(reader->GetOutputPort());
    flip->SetResliceAxesDirectionCosines(
      -1,0,0, 0,-1,0, 0,0,1);
    flip->SetOutputSpacing(spacing);
    flip->Update();

    image = flip->GetOutput();
    }

  // get the data
  data->CopyStructure(image);
  data->GetPointData()->PassData(image->GetPointData());

  if (coordSystem == DICOMCoords)
    {
    // generate the matrix, but modify to use DICOM coords
    static double xyFlipMatrix[16] =
      { -1, 0, 0, 0,  0, -1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1 };
    // correct for the flip that was done earlier
    vtkMatrix4x4::Multiply4x4(*reader->GetDirectionCosines()->Element,
                              xyFlipMatrix, *matrix->Element);
    // do the left/right, up/down dicom-to-minc transformation
    vtkMatrix4x4::Multiply4x4(xyFlipMatrix, *matrix->Element, *matrix->Element);
    matrix->Modified();
    }
  else
    {
    matrix->DeepCopy(reader->GetDirectionCosines());
    }
}

#ifdef AIRS_USE_NIFTI
void ReadNIFTIImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *fileName,
  int coordSystem)
{
  // read the image
  vtkSmartPointer<vtkNIFTIReader> reader =
    vtkSmartPointer<vtkNIFTIReader>::New();

  reader->SetFileName(fileName);
  reader->Update();

  vtkSmartPointer<vtkImageData> image = reader->GetOutput();

  if (coordSystem == DICOMCoords)
    {
    double spacing[3];
    reader->GetOutput()->GetSpacing(spacing);
    spacing[0] = fabs(spacing[0]);
    spacing[1] = fabs(spacing[1]);
    spacing[2] = fabs(spacing[2]);

    // flip the image rows into a DICOM-style ordering
    vtkSmartPointer<vtkImageReslice> flip =
      vtkSmartPointer<vtkImageReslice>::New();

    flip->SetInputConnection(reader->GetOutputPort());
    flip->SetResliceAxesDirectionCosines(
      -1,0,0, 0,-1,0, 0,0,1);
    flip->SetOutputSpacing(spacing);
    flip->Update();

    image = flip->GetOutput();
    }

  // get the data
  data->CopyStructure(image);
  data->GetPointData()->PassData(image->GetPointData());

  // get the SForm or QForm matrix if present
  static double nMatrix[16] =
    { 1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1 };
  if (reader->GetSFormMatrix())
    {
    vtkMatrix4x4::DeepCopy(nMatrix, reader->GetSFormMatrix());
    }
  else if (reader->GetQFormMatrix())
    {
    vtkMatrix4x4::DeepCopy(nMatrix, reader->GetQFormMatrix());
    }

  if (coordSystem == DICOMCoords)
    {
    // generate the matrix, but modify to use DICOM coords
    static double xyFlipMatrix[16] =
      { -1, 0, 0, 0,  0, -1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1 };
    // correct for the flip that was done earlier
    vtkMatrix4x4::Multiply4x4(nMatrix, xyFlipMatrix, *matrix->Element);
    // do the left/right, up/down dicom-to-minc transformation
    vtkMatrix4x4::Multiply4x4(xyFlipMatrix, *matrix->Element, *matrix->Element);
    matrix->Modified();
    }
  else
    {
    matrix->DeepCopy(nMatrix);
    }
}
#endif /* AIRS_USE_NIFTI */

void ReadImage(
  vtkImageData *image, vtkMatrix4x4 *matrix,
  const char *filename, int coordSystem)
{
  size_t n = strlen(filename);
  if (n > 4 && strcmp(&filename[n-4], ".mnc") == 0)
    {
    ReadMINCImage(image, matrix, filename, coordSystem);
    }
  else if ((n > 4 && strcmp(&filename[n-4], ".nii") == 0) ||
           (n > 7 && strcmp(&filename[n-7], ".nii.gz") == 0))
    {
#ifdef AIRS_USE_NIFTI
    ReadNIFTIImage(image, matrix, filename, coordSystem);
#else
    fprintf(stderr, "NIFTI files are not supported.\n");
    exit(1);
#endif
    }
  else
    {
    ReadDICOMImage(image, matrix, filename, coordSystem);
    }
}

int CoordSystem(const char *filename)
{
  size_t n = strlen(filename);
  if ((n > 4 && strcmp(&filename[n-4], ".mnc") == 0) ||
      (n > 4 && strcmp(&filename[n-4], ".nii") == 0) ||
      (n > 7 && strcmp(&filename[n-7], ".nii.gz") == 0))
    {
    return NIFTICoords;
    }

  return DICOMCoords;
}

void SetViewFromMatrix(
  vtkRenderer *renderer,
  vtkInteractorStyleImage *istyle,
  vtkMatrix4x4 *matrix,
  int coordSystem)
{
  istyle->SetCurrentRenderer(renderer);

  // This view assumes the data uses the DICOM Patient Coordinate System.
  // It provides a right-is-left view of axial and coronal images
  double viewRight[4] = { 1.0, 0.0, 0.0, 0.0 };
  double viewUp[4] = { 0.0, 1.0, 0.0, 0.0 };

  if (coordSystem == DICOMCoords)
    {
    viewUp[1] = -1.0;
    }

  matrix->MultiplyPoint(viewRight, viewRight);
  matrix->MultiplyPoint(viewUp, viewUp);

  istyle->SetImageOrientation(viewRight, viewUp);
}

void ReadMatrix(vtkMatrix4x4 *matrix, const char *xfminput)
{
  size_t l = strlen(xfminput);
  if (l >= 4 && strcmp(xfminput + (l - 4), ".xfm") == 0)
    {
    // MNI transform file (always in RAS coords)
    vtkSmartPointer<vtkMNITransformReader> reader =
      vtkSmartPointer<vtkMNITransformReader>::New();
    reader->SetFileName(xfminput);
    reader->Update();
    vtkLinearTransform *transform =
      vtkLinearTransform::SafeDownCast(reader->GetTransform());
    if (transform)
      {
      matrix->DeepCopy(transform->GetMatrix());
      }
    }
  else if ((l >= 4 && strcmp(xfminput + (l - 4), ".txt") == 0) ||
           (l >= 4 && strcmp(xfminput + (l - 4), ".tfm") == 0))
    {
    // ITK transform file (always in DICOM coords)
    vtkSmartPointer<vtkITKXFMReader> reader =
      vtkSmartPointer<vtkITKXFMReader>::New();
    reader->SetFileName(xfminput);
    reader->Update();
    vtkLinearTransform *transform =
      vtkLinearTransform::SafeDownCast(reader->GetTransform());
    if (transform)
      {
      matrix->DeepCopy(transform->GetMatrix());
      }
    }
  else
    {
    // Space-delimited text file
    double elements[16] = {
      1.0, 0.0, 0.0, 0.0,
      0.0, 1.0, 0.0, 0.0,
      0.0, 0.0, 1.0, 0.0,
      0.0, 0.0, 0.0, 1.0 };

    ifstream infile(xfminput);
    for (int i = 0; infile && i < 16; i++)
      {
      infile >> elements[i];
      }
    infile.close();
    matrix->DeepCopy(elements);
    }
}

void WriteMatrix(
  vtkMatrix4x4 *matrix, const char *xfmfile, const double center[3])
{
  vtkSmartPointer<vtkTransform> transform =
    vtkSmartPointer<vtkTransform>::New();
  transform->Concatenate(matrix);

  size_t l = strlen(xfmfile);
  if (l >= 4 && strcmp(xfmfile + (l - 4), ".xfm") == 0)
    {
    // MNI transform file (always in RAS coords)
    vtkSmartPointer<vtkMNITransformWriter> writer =
      vtkSmartPointer<vtkMNITransformWriter>::New();
    writer->SetFileName(xfmfile);
    writer->SetTransform(transform);
    writer->Update();
    }
  else if ((l >= 4 && strcmp(xfmfile + (l - 4), ".txt") == 0) ||
           (l >= 4 && strcmp(xfmfile + (l - 4), ".tfm") == 0))
    {
    // ITK transform file (always in DICOM coords)
    vtkSmartPointer<vtkITKXFMWriter> writer =
      vtkSmartPointer<vtkITKXFMWriter>::New();
    writer->SetFileName(xfmfile);
    writer->SetTransform(transform);
    writer->SetTransformCenter(center);
    writer->Write();
    }
  else
    {
    // Delimited text file
    ofstream outfile(xfmfile, ios::out);
    for (int i = 0; i < 4; i++)
      {
      outfile << matrix->Element[i][0] << "  "
              << matrix->Element[i][1] << "  "
              << matrix->Element[i][2] << "  "
              << matrix->Element[i][3] << "\n";
      }
    outfile.close();
    }
}

void WriteScreenshot(vtkWindow *window, const char *filename)
{
  vtkSmartPointer<vtkWindowToImageFilter> snap =
    vtkSmartPointer<vtkWindowToImageFilter>::New();
  snap->SetInput(window);
  snap->Update();

  size_t l = strlen(filename);
  if (l >= 4 && strcmp(filename + (l - 4), ".png") == 0)
    {
    vtkSmartPointer<vtkPNGWriter> snapWriter =
      vtkSmartPointer<vtkPNGWriter>::New();
    snapWriter->SetInputConnection(snap->GetOutputPort());
    snapWriter->SetFileName(filename);
    snapWriter->Write();
    }
  else if ((l >= 4 && strcmp(filename + (l - 4), ".jpg") == 0) ||
           (l >= 5 && strcmp(filename + (l - 5), ".jpeg") == 0))
    {
    vtkSmartPointer<vtkJPEGWriter> snapWriter =
      vtkSmartPointer<vtkJPEGWriter>::New();
    snapWriter->SetInputConnection(snap->GetOutputPort());
    snapWriter->SetFileName(filename);
    snapWriter->Write();
    }
  else if ((l >= 4 && strcmp(filename + (l - 4), ".tif") == 0) ||
           (l >= 5 && strcmp(filename + (l - 5), ".tiff") == 0))
    {
    vtkSmartPointer<vtkTIFFWriter> snapWriter =
      vtkSmartPointer<vtkTIFFWriter>::New();
    snapWriter->SetInputConnection(snap->GetOutputPort());
    snapWriter->SetFileName(filename);
    snapWriter->Write();
    }
}

};

struct register_options
{
  int dimensionality;  // -D --dimensionality
  int metric;          // -M --metric
  int transform;       // -T --transform
  int coords;          // -C --coords
  int interactive;     // -I --interactive
  int checkonly;       // -J --checkonly
  const char *initial; // -i --initial (initial transform)
  const char *output;  // -o (output transform)
  const char *screenshot; // -s (output screenshot)
  const char *source;
  const char *target;
};

void register_initialize_options(register_options *options)
{
  options->dimensionality = 3;
  options->metric = vtkImageRegistration::NormalizedMutualInformation;
  options->transform = vtkImageRegistration::Rigid;
  options->coords = NativeCoords;
  options->checkonly = 0;
  options->interactive = 0;
  options->initial = NULL;
  options->screenshot = NULL;
  options->output = NULL;
  options->source = NULL;
  options->target = NULL;
}

const char *check_next_arg(
  int argc, char *argv[], int *argi, const char *possib[])
{
  const char *op = argv[*argi - 1];
  if (*argi >= argc ||
      argv[*argi][0] == '-')
    {
    fprintf(stderr, "The option \"%s\" must be followed by an argument\n", op);
    exit(1);
    }
  const char *arg = argv[(*argi)++];

  if (possib == 0)
    {
    return arg;
    }

  for (const char **t = possib; *t != 0; t++)
    {
    if (strcmp(*t, arg) == 0)
      {
      return arg;
      }
    }

  fprintf(stderr, "Incorrect value for option \"%s\": %s\n",
          op, arg);
  fprintf(stderr, "Allowed values:");
  for (const char **u = possib; *u != 0; u++)
    {
    fprintf(stderr, "%s", *u);
    }
  fprintf(stderr, "\n");
  exit(1);

  return 0;
}

int register_read_options(
  int argc, char *argv[], register_options *options)
{
  static const char *dimensionality_args[] = {
    "2", "3", 0 };
  static const char *metric_args[] = {
    "SquaredDifference", "SD",
    "CrossCorrelation", "CC",
    "NormalizedCrossCorrelation", "NCC",
    "NeighborhoodCorrelation", "NC",
    "MutualInformation", "MI",
    "NormalizedMutualInformation", "NMI",
    0 };
  static const char *transform_args[] = {
    "Translation", "TR",
    "Rigid", "RI",
    "Similarity", "SI",
    "ScaleSourceAxes", "SS",
    "ScaleTargetAxes", "ST",
    "Affine", "AF",
    0 };
  static const char *coords_args[] = {
    "DICOM", "LPS",
    "NIFTI", "MINC", "RAS",
    0 };

  int argi = 1;
  while (argi < argc)
    {
    const char *arg = argv[argi++];
    if (arg[0] != '-')
      {
      if (options->source == 0)
        {
        options->source = arg;
        }
      else if (options->target == 0)
        {
        options->target = arg;
        }
      else
        {
        fprintf(stderr, "Too many files listed on command line\n");
        exit(1);
        }
      }
    else
      {
      if (strcmp(arg, "-D") == 0 ||
          strcmp(arg, "--dimensionality") == 0)
        {
        arg = check_next_arg(argc, argv, &argi, dimensionality_args);
        options->dimensionality = (arg[0] == '2' ? 2 : 3);
        }
      else if (strcmp(arg, "-M") == 0 ||
               strcmp(arg, "--metric") == 0)
        {
        arg = check_next_arg(argc, argv, &argi, metric_args);
        if (strcmp(arg, "SquaredDifference") == 0 ||
            strcmp(arg, "SD") == 0)
          {
          options->metric = vtkImageRegistration::SquaredDifference;
          }
        else if (strcmp(arg, "CrossCorrelation") == 0 ||
                 strcmp(arg, "CC") == 0)
          {
          options->metric = vtkImageRegistration::CrossCorrelation;
          }
        else if (strcmp(arg, "NormalizedCrossCorrelation") == 0 ||
                 strcmp(arg, "NCC") == 0)
          {
          options->metric = vtkImageRegistration::NormalizedCrossCorrelation;
          }
        else if (strcmp(arg, "NeighborhoodCorrelation") == 0 ||
                 strcmp(arg, "NC") == 0)
          {
          options->metric = vtkImageRegistration::NeighborhoodCorrelation;
          }
        else if (strcmp(arg, "MutualInformation") == 0 ||
                 strcmp(arg, "MI") == 0)
          {
          options->metric = vtkImageRegistration::MutualInformation;
          }
        else if (strcmp(arg, "NormalizedMutualInformation") == 0 ||
                 strcmp(arg, "NMI") == 0)
          {
          options->metric = vtkImageRegistration::NormalizedMutualInformation;
          }
        }
      else if (strcmp(arg, "-T") == 0 ||
               strcmp(arg, "--transform") == 0)
        {
        arg = check_next_arg(argc, argv, &argi, transform_args);
        if (strcmp(arg, "Translation") == 0 ||
            strcmp(arg, "TR") == 0)
          {
          options->transform = vtkImageRegistration::Translation;
          }
        else if (strcmp(arg, "Rigid") == 0 ||
                 strcmp(arg, "RI") == 0)
          {
          options->transform = vtkImageRegistration::Rigid;
          }
        else if (strcmp(arg, "Similarity") == 0 ||
                 strcmp(arg, "SI") == 0)
          {
          options->transform = vtkImageRegistration::Similarity;
          }
        else if (strcmp(arg, "ScaleSourceAxes") == 0 ||
                 strcmp(arg, "SS") == 0)
          {
          options->transform = vtkImageRegistration::ScaleSourceAxes;
          }
        else if (strcmp(arg, "ScaleTargetAxes") == 0 ||
                 strcmp(arg, "ST") == 0)
          {
          options->transform = vtkImageRegistration::ScaleTargetAxes;
          }
        else if (strcmp(arg, "Affine") == 0 ||
                 strcmp(arg, "AF") == 0)
          {
          options->transform = vtkImageRegistration::Affine;
          }
        }
      else if (strcmp(arg, "-C") == 0 ||
               strcmp(arg, "--coords") == 0)
        {
        arg = check_next_arg(argc, argv, &argi, coords_args);
        if (strcmp(arg, "DICOM") == 0 ||
            strcmp(arg, "LPS") == 0)
          {
          options->coords = DICOMCoords;
          }
        else if (strcmp(arg, "MINC") == 0 ||
                 strcmp(arg, "NIFTI") == 0 ||
                 strcmp(arg, "RAS") == 0)
          {
          options->coords = NIFTICoords;
          }
        }
      else if (strcmp(arg, "-I") == 0 ||
               strcmp(arg, "--interactive") == 0)
        {
        options->interactive = 1;
        }
      else if (strcmp(arg, "-J") == 0 ||
               strcmp(arg, "--checkonly") == 0)
        {
        options->interactive = 1;
        options->checkonly = 1;
        }
      else if (strcmp(arg, "-i") == 0 ||
               strcmp(arg, "--initial") == 0)
        {
        arg = check_next_arg(argc, argv, &argi, 0);
        options->initial = arg;
        }
      else if (strcmp(arg, "-s") == 0 ||
               strcmp(arg, "--screenshot") == 0)
        {
        arg = check_next_arg(argc, argv, &argi, 0);
        options->screenshot = arg;
        }
      else if (strcmp(arg, "-o") == 0)
        {
        arg = check_next_arg(argc, argv, &argi, 0);
        options->output = arg;
        }
      else
        {
        fprintf(stderr, "Unrecognized option \"%s\"\n", arg);
        exit(1);
        }
      }
    }

  return 1;
}

int main(int argc, char *argv[])
{
  register_options options;
  register_initialize_options(&options);
  register_read_options(argc, argv, &options);

  // -------------------------------------------------------
  // the files
  const char *xfminput = options.initial;
  const char *xfmfile = options.output;
  const char *sourcefile = options.source;
  const char *targetfile = options.target;
  bool display = (options.interactive != 0 ||
                  options.screenshot != 0);

  const char *xfiles[2];
  xfiles[0] = xfminput;
  xfiles[1] = xfmfile;
  for (int xi = 0; xi < 2; xi++)
    {
    const char *xfile = xfiles[xi];
    if (xfile != 0)
      {
      size_t m = strlen(xfile);
      if (m < 4 ||
          (strcmp(&xfile[m-4], ".xfm") != 0) &&
          (strcmp(&xfile[m-4], ".tfm") != 0) &&
          (strcmp(&xfile[m-4], ".txt") != 0) &&
          (strcmp(&xfile[m-4], ".mat") != 0))
        {
        fprintf(stderr, "Transform file must end in .xfm, .tfm, or .txt\n");
        return 1;
        }
      }
    }

  // -------------------------------------------------------
  // parameters for registration

  int interpolatorType = vtkImageRegistration::Linear;
  double transformTolerance = 0.1; // tolerance on transformation result
  int numberOfBins = 64; // for Mattes' mutual information
  double initialBlurFactor = 4.0;

  // -------------------------------------------------------
  // load the images

  if (options.coords == NativeCoords)
    {
    int ic = CoordSystem(sourcefile);
    int oc = CoordSystem(targetfile);

    if (ic == DICOMCoords || oc == DICOMCoords)
      {
      options.coords = DICOMCoords;
      }
    else
      {
      options.coords = NIFTICoords;
      }
    }

  vtkSmartPointer<vtkImageData> sourceImage =
    vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkMatrix4x4> sourceMatrix =
    vtkSmartPointer<vtkMatrix4x4>::New();
  ReadImage(sourceImage, sourceMatrix, sourcefile, options.coords);

  vtkSmartPointer<vtkImageData> targetImage =
    vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkMatrix4x4> targetMatrix =
    vtkSmartPointer<vtkMatrix4x4>::New();
  ReadImage(targetImage, targetMatrix, targetfile, options.coords);

  // -------------------------------------------------------
  // save the original source matrix
  vtkSmartPointer<vtkMatrix4x4> originalSourceMatrix =
    vtkSmartPointer<vtkMatrix4x4>::New();
  originalSourceMatrix->DeepCopy(sourceMatrix);

  // -------------------------------------------------------
  // load the initial matrix
  if (xfminput)
    {
    vtkSmartPointer<vtkMatrix4x4> initialMatrix =
      vtkSmartPointer<vtkMatrix4x4>::New();
    ReadMatrix(initialMatrix, xfminput);
    vtkMatrix4x4::Multiply4x4(initialMatrix, sourceMatrix, sourceMatrix);
    }

  // -------------------------------------------------------
  // display the images

  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindowInteractor> interactor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  vtkSmartPointer<vtkInteractorStyleImage> istyle =
    vtkSmartPointer<vtkInteractorStyleImage>::New();

  istyle->SetInteractionModeToImageSlicing();
  interactor->SetInteractorStyle(istyle);
  renderWindow->SetInteractor(interactor);
  renderWindow->AddRenderer(renderer);

  vtkSmartPointer<vtkImageSlice> sourceActor =
    vtkSmartPointer<vtkImageSlice>::New();
  vtkSmartPointer<vtkImageResliceMapper> sourceMapper =
    vtkSmartPointer<vtkImageResliceMapper>::New();
  vtkSmartPointer<vtkImageProperty> sourceProperty =
    vtkSmartPointer<vtkImageProperty>::New();

  sourceMapper->SetInput(sourceImage);
  sourceMapper->SliceAtFocalPointOn();
  sourceMapper->SliceFacesCameraOn();
  sourceMapper->ResampleToScreenPixelsOff();

  double sourceRange[2];
  sourceImage->GetScalarRange(sourceRange);
  sourceProperty->SetInterpolationTypeToLinear();
  sourceProperty->SetColorWindow((sourceRange[1]-sourceRange[0]));
  sourceProperty->SetColorLevel(0.5*(sourceRange[0]+sourceRange[1]));
  sourceProperty->CheckerboardOn();

  sourceActor->SetMapper(sourceMapper);
  sourceActor->SetProperty(sourceProperty);
  sourceActor->SetUserMatrix(sourceMatrix);

  vtkSmartPointer<vtkImageSlice> targetActor =
    vtkSmartPointer<vtkImageSlice>::New();
  vtkSmartPointer<vtkImageResliceMapper> targetMapper =
    vtkSmartPointer<vtkImageResliceMapper>::New();
  vtkSmartPointer<vtkImageProperty> targetProperty =
    vtkSmartPointer<vtkImageProperty>::New();

  targetMapper->SetInput(targetImage);
  targetMapper->SliceAtFocalPointOn();
  targetMapper->SliceFacesCameraOn();
  targetMapper->ResampleToScreenPixelsOff();

  double targetRange[2];
  targetImage->GetScalarRange(targetRange);
  targetProperty->SetInterpolationTypeToLinear();
  targetProperty->SetColorWindow((targetRange[1]-targetRange[0]));
  targetProperty->SetColorLevel(0.5*(targetRange[0]+targetRange[1]));

  targetActor->SetMapper(targetMapper);
  targetActor->SetProperty(targetProperty);
  targetActor->SetUserMatrix(targetMatrix);

  vtkSmartPointer<vtkImageStack> imageStack =
    vtkSmartPointer<vtkImageStack>::New();
  imageStack->AddImage(targetActor);
  imageStack->AddImage(sourceActor);

  renderer->AddViewProp(imageStack);
  renderer->SetBackground(0,0,0);

  renderWindow->SetSize(1024,1024);

  double bounds[6], center[4], tspacing[3];
  int extent[6];
  targetImage->GetBounds(bounds);
  targetImage->GetExtent(extent);
  targetImage->GetSpacing(tspacing);
  center[0] = 0.5*(bounds[0] + bounds[1]);
  center[1] = 0.5*(bounds[2] + bounds[3]);
  center[2] = 0.5*(bounds[4] + bounds[5]);
  center[3] = 1.0;
  targetMatrix->MultiplyPoint(center, center);

  vtkCamera *camera = renderer->GetActiveCamera();
  renderer->ResetCamera();
  camera->SetFocalPoint(center);
  camera->ParallelProjectionOn();
  camera->SetParallelScale(0.5*(bounds[3] - bounds[2]));
  SetViewFromMatrix(renderer, istyle, targetMatrix, options.coords);
  renderer->ResetCameraClippingRange();

  double checkSpacing = (extent[3] - extent[2] + 7)/7*tspacing[1];
  sourceProperty->SetCheckerboardSpacing(checkSpacing, checkSpacing);

  if (display)
    {
    renderWindow->Render();
    }

  // -------------------------------------------------------
  // prepare for registration

  // get information about the images
  double targetSpacing[3], sourceSpacing[3];
  targetImage->GetSpacing(targetSpacing);
  sourceImage->GetSpacing(sourceSpacing);

  for (int jj = 0; jj < 3; jj++)
    {
    targetSpacing[jj] = fabs(targetSpacing[jj]);
    sourceSpacing[jj] = fabs(sourceSpacing[jj]);
    }

  double minSpacing = sourceSpacing[0];
  if (minSpacing > sourceSpacing[1])
    {
    minSpacing = sourceSpacing[1];
    }
  if (minSpacing > sourceSpacing[2])
    {
    minSpacing = sourceSpacing[2];
    }

  // blur source image with Blackman-windowed sinc
  vtkSmartPointer<vtkImageSincInterpolator> sourceBlurKernel =
    vtkSmartPointer<vtkImageSincInterpolator>::New();
  sourceBlurKernel->SetWindowFunctionToBlackman();

  // reduce the source resolution
  vtkSmartPointer<vtkImageResize> sourceBlur =
    vtkSmartPointer<vtkImageResize>::New();
  sourceBlur->SetInput(sourceImage);
  sourceBlur->SetResizeMethodToOutputSpacing();
  sourceBlur->SetInterpolator(sourceBlurKernel);
  sourceBlur->InterpolateOn();

  // blur target with Blackman-windowed sinc
  vtkSmartPointer<vtkImageSincInterpolator> targetBlurKernel =
    vtkSmartPointer<vtkImageSincInterpolator>::New();
  targetBlurKernel->SetWindowFunctionToBlackman();

  // keep target at full resolution
  vtkSmartPointer<vtkImageResize> targetBlur =
    vtkSmartPointer<vtkImageResize>::New();
  targetBlur->SetInput(targetImage);
  targetBlur->SetResizeMethodToOutputSpacing();
  targetBlur->SetInterpolator(targetBlurKernel);
  targetBlur->InterpolateOn();

  // get the initial transformation
  vtkSmartPointer<vtkMatrix4x4> matrix =
    vtkSmartPointer<vtkMatrix4x4>::New();
  matrix->DeepCopy(targetMatrix);
  matrix->Invert();
  vtkMatrix4x4::Multiply4x4(matrix, sourceMatrix, matrix);

  // set up the registration
  vtkSmartPointer<vtkImageRegistration> registration =
    vtkSmartPointer<vtkImageRegistration>::New();
  registration->SetTargetImageInputConnection(targetBlur->GetOutputPort());
  registration->SetSourceImageInputConnection(sourceBlur->GetOutputPort());
  registration->SetTransformDimensionality(options.dimensionality);
  registration->SetTransformType(options.transform);
  registration->SetMetricType(options.metric);
  registration->SetInterpolatorType(interpolatorType);
  registration->SetJointHistogramSize(numberOfBins,numberOfBins);
  registration->SetMetricTolerance(1e-4);
  registration->SetTransformTolerance(transformTolerance);
  registration->SetMaximumNumberOfIterations(500);
  if (xfminput)
    {
    registration->SetInitializerTypeToNone();
    }
  else
    {
    registration->SetInitializerTypeToCentered();
    }
  registration->Initialize(matrix);

  // -------------------------------------------------------
  // make a timer
  vtkSmartPointer<vtkTimerLog> timer =
    vtkSmartPointer<vtkTimerLog>::New();
  double startTime = timer->GetUniversalTime();
  double lastTime = startTime;

  // -------------------------------------------------------
  // do the registration

  // the registration starts at low-resolution
  double blurFactor = initialBlurFactor;
  // two stages for each resolution:
  // first without interpolation, and then with interpolation
  int stage = 0;
  // will be set to "true" when registration is initialized
  bool initialized = false;

  while (options.checkonly == 0)
    {
    if (stage == 0)
      {
      registration->SetInterpolatorTypeToNearest();
      registration->SetTransformTolerance(minSpacing*blurFactor);
      }
    else
      {
      registration->SetInterpolatorType(interpolatorType);
      registration->SetTransformTolerance(transformTolerance*blurFactor);
      }
    if (blurFactor < 1.1)
      {
      // full resolution: no blurring or resampling
      sourceBlur->SetInterpolator(0);
      sourceBlur->InterpolateOff();
      sourceBlur->SetOutputSpacing(sourceSpacing);
      sourceBlur->Update();

      targetBlur->SetInterpolator(0);
      sourceBlur->InterpolateOff();
      targetBlur->SetOutputSpacing(targetSpacing);
      targetBlur->Update();
      }
    else
      {
      // reduced resolution: set the blurring
      double spacing[3];
      for (int j = 0; j < 3; j++)
        {
        spacing[j] = blurFactor*minSpacing;
        if (spacing[j] < sourceSpacing[j])
          {
          spacing[j] = sourceSpacing[j];
          }
        }

      sourceBlurKernel->SetBlurFactors(
        spacing[0]/sourceSpacing[0],
        spacing[1]/sourceSpacing[1],
        spacing[2]/sourceSpacing[2]);

      sourceBlur->SetOutputSpacing(spacing);
      sourceBlur->Update();

      targetBlurKernel->SetBlurFactors(
        blurFactor*minSpacing/targetSpacing[0],
        blurFactor*minSpacing/targetSpacing[1],
        blurFactor*minSpacing/targetSpacing[2]);

      targetBlur->Update();
      }

    if (initialized)
      {
      // re-initialize with the matrix from the previous step
      registration->SetInitializerTypeToNone();
      matrix->DeepCopy(registration->GetTransform()->GetMatrix());
      }

    registration->Initialize(matrix);

    initialized = true;

    while (registration->Iterate())
      {
      // registration->UpdateRegistration();
      // will iterate until convergence or failure
      vtkMatrix4x4::Multiply4x4(
        targetMatrix,registration->GetTransform()->GetMatrix(),sourceMatrix);
      sourceMatrix->Modified();
      if (display)
        {
        interactor->Render();
        }
      }

    double newTime = timer->GetUniversalTime();
    cout << "blur " << blurFactor << " stage " << stage << " took "
         << (newTime - lastTime) << "s and "
         << registration->GetNumberOfEvaluations() << " evaluations" << endl;
    lastTime = newTime;

    // prepare for next iteration
    if (stage == 1)
      {
      blurFactor /= 2.0;
      if (blurFactor < 0.9)
        {
        break;
        }
      }
    stage = (stage + 1) % 2;
    }

  cout << "registration took " << (lastTime - startTime) << "s" << endl;

  // -------------------------------------------------------
  // write the output matrix
  if (xfmfile)
    {
    vtkMatrix4x4 *rmatrix = registration->GetTransform()->GetMatrix();
    vtkSmartPointer<vtkMatrix4x4> wmatrix =
      vtkSmartPointer<vtkMatrix4x4>::New();
    wmatrix->DeepCopy(originalSourceMatrix);
    wmatrix->Invert();
    vtkMatrix4x4::Multiply4x4(rmatrix, wmatrix, wmatrix);
    vtkMatrix4x4::Multiply4x4(targetMatrix, wmatrix, wmatrix);

    WriteMatrix(wmatrix, xfmfile, center);
    }

  // -------------------------------------------------------
  // capture a screen shot
  if (options.screenshot)
    {
    WriteScreenshot(renderWindow, options.screenshot);
    }

  // -------------------------------------------------------
  // allow user to interact

  if (options.interactive)
    {
    interactor->Start();
    }

  return 0;
}
