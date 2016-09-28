/*=========================================================================

Program:   Atamai Image Registration and Segmentation
Module:    skullstrip.cxx

   This software is distributed WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/

// Apply Stephen Smith's BET algorithm to an MRI image.

#include <vtkSmartPointer.h>

#include <vtkImageReslice.h>
#include <vtkImageHistogramStatistics.h>
#include <vtkROIStencilSource.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkLookupTable.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>
#include <vtkMath.h>
#include <vtkCommand.h>

#include <vtkMINCImageReader.h>
#include <vtkMINCImageWriter.h>
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
#include <vtkSTLWriter.h>
#include <vtkMNIObjectWriter.h>
#include <vtkPolyDataWriter.h>

#include <vtkTimerLog.h>
#include <vtkVersion.h>

#include <vtksys/SystemTools.hxx>

#include "AIRSConfig.h"
#include "vtkImageMRIBrainExtractor.h"

// optional readers
#ifdef AIRS_USE_DICOM
#define AIRS_USE_NIFTI
#include <vtkNIFTIReader.h>
#include <vtkNIFTIWriter.h>
#include <vtkDICOMReader.h>
#include <vtkDICOMSorter.h>
#include <vtkDICOMMRGenerator.h>
#include <vtkDICOMCTGenerator.h>
#include <vtkDICOMWriter.h>
#include <vtkDICOMMetaData.h>
#include <vtkGlobFileNames.h>
#endif

#include <vector>
#include <string>

#include <stdlib.h>

// A macro to assist VTK 5 backwards compatibility
#if VTK_MAJOR_VERSION >= 6
#define SET_INPUT_DATA SetInputData
#define SET_STENCIL_DATA SetStencilData
#else
#define SET_INPUT_DATA SetInput
#define SET_STENCIL_DATA SetStencil
#endif

// coord systems
enum { NativeCoords, DICOMCoords, NIFTICoords };

// file types
enum { DICOMImage, NIFTIImage, MINCImage, LastImageType = MINCImage,
       STLSurface, OBJSurface, VTKSurface, LastSurfaceType = VTKSurface };

// internal methods for reading images, these methods read the image
// into the specified data object and also provide a matrix for converting
// the data coordinates into patient coordinates.
namespace {

int GuessFileType(const char *filename)
{
  size_t n = strlen(filename);

  if (n > 4 && strcmp(&filename[n-4], ".mnc") == 0)
    {
    return MINCImage;
    }
  if ((n > 4 && strcmp(&filename[n-4], ".nii") == 0) ||
      (n > 7 && strcmp(&filename[n-7], ".nii.gz") == 0))
    {
    return NIFTIImage;
    }
  if (n > 4 && strcmp(&filename[n-4], ".stl") == 0)
    {
    return STLSurface;
    }
  if (n > 4 && strcmp(&filename[n-4], ".obj") == 0)
    {
    return OBJSurface;
    }
  if (n > 4 && strcmp(&filename[n-4], ".vtk") == 0)
    {
    return VTKSurface;
    }

  return DICOMImage;
}

#ifdef AIRS_USE_DICOM
vtkDICOMReader *ReadDICOMImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *directoryName,
  int coordSystem)
{
  // get the files
  std::string dirString = directoryName;
  vtksys::SystemTools::ConvertToUnixSlashes(dirString);
  vtkSmartPointer<vtkGlobFileNames> glob =
    vtkSmartPointer<vtkGlobFileNames>::New();
  glob->SetDirectory(dirString.c_str());
  glob->AddFileNames("*");

  // sort the files
  vtkSmartPointer<vtkDICOMSorter> sorter =
    vtkSmartPointer<vtkDICOMSorter>::New();
  sorter->SetInputFileNames(glob->GetFileNames());
  sorter->Update();

  if (sorter->GetNumberOfSeries() == 0)
    {
    fprintf(stderr, "Folder contains no DICOM files: %s\n", directoryName);
    exit(1);
    }
  else if (sorter->GetNumberOfSeries() > 1)
    {
    fprintf(stderr, "Folder contains more than one DICOM series: %s\n",
            directoryName);
    exit(1);
    }

  // read the image
  vtkDICOMReader *reader = vtkDICOMReader::New();
  reader->SetFileNames(sorter->GetFileNamesForSeries(0));

  if (coordSystem == NIFTICoords)
    {
    reader->SetMemoryRowOrderToBottomUp();
    }
  else
    {
    reader->SetMemoryRowOrderToFileNative();
    }

  reader->UpdateInformation();
  if (reader->GetErrorCode())
    {
    exit(1);
    }

  // when reading images, only read 1st component if the
  // image has multiple components or multiple time points
  vtkIntArray *fileArray = reader->GetFileIndexArray();

  // create a filtered list of files
  vtkSmartPointer<vtkStringArray> fileNames =
    vtkSmartPointer<vtkStringArray>::New();
  vtkIdType n = fileArray->GetNumberOfTuples();
  for (vtkIdType i = 0; i < n; i++)
    {
    fileNames->InsertNextValue(
      reader->GetFileNames()->GetValue(fileArray->GetComponent(i, 0)));
    }
  reader->SetDesiredTimeIndex(0);
  reader->SetFileNames(fileNames);

  reader->Update();
  if (reader->GetErrorCode())
    {
    exit(1);
    }

  vtkImageData *image = reader->GetOutput();

  // get the data
  data->CopyStructure(image);
  data->GetPointData()->PassData(image->GetPointData());

  // get the matrix
  matrix->DeepCopy(reader->GetPatientMatrix());

  return reader;
}

void WriteDICOMImage(
  vtkImageReader2 *sourceReader, vtkImageReader2 *targetReader,
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *directoryName,
  int vtkNotUsed(coordSystem))
{
  if (vtksys::SystemTools::FileExists(directoryName))
    {
    if (!vtksys::SystemTools::FileIsDirectory(directoryName))
      {
      fprintf(stderr, "option -o must give a DICOM directory, not a file.\n");
      exit(1);
      }
    }
  else if (!vtksys::SystemTools::MakeDirectory(directoryName))
    {
    fprintf(stderr, "Cannot create directory: %s\n", directoryName);
    exit(1);
    }

  // get the meta data
  vtkDICOMReader *reader = vtkDICOMReader::SafeDownCast(targetReader);
  vtkDICOMReader *reader2 = vtkDICOMReader::SafeDownCast(sourceReader);

  vtkSmartPointer<vtkDICOMMetaData> meta =
    vtkSmartPointer<vtkDICOMMetaData>::New();

  if (reader)
    {
    // copy the bulk of the meta data from the target image
    meta->DeepCopy(reader->GetMetaData());
    meta->SetAttributeValue(DC::SeriesNumber,
      meta->GetAttributeValue(DC::SeriesNumber).AsUnsignedInt() + 1000);
    std::string seriesDescription =
      meta->GetAttributeValue(DC::SeriesDescription).AsString() + " SEG";
    if (seriesDescription.size() < 64)
      {
      meta->SetAttributeValue(DC::SeriesDescription, seriesDescription);
      }
    }
  if (reader2)
    {
    // set the frame of reference from the source image
    meta->SetAttributeValue(DC::FrameOfReferenceUID,
      reader2->GetMetaData()->GetAttributeValue(
      DC::FrameOfReferenceUID));
    }

  // make the generator
  vtkSmartPointer<vtkDICOMMRGenerator> mrgenerator =
    vtkSmartPointer<vtkDICOMMRGenerator>::New();
  vtkSmartPointer<vtkDICOMCTGenerator> ctgenerator =
    vtkSmartPointer<vtkDICOMCTGenerator>::New();
  vtkDICOMGenerator *generator = 0;
  if (reader)
    {
    std::string SOPClass =
      meta->GetAttributeValue(DC::SOPClassUID).AsString();
    if (SOPClass == "1.2.840.10008.5.1.4.1.1.2")
      {
      generator = ctgenerator;
      }
    else if (SOPClass == "1.2.840.10008.5.1.4.1.1.4")
      {
      generator = mrgenerator;
      }
    }

  // prepare the writer to write the image
  vtkSmartPointer<vtkDICOMWriter> writer =
    vtkSmartPointer<vtkDICOMWriter>::New();
  if (generator)
    {
    writer->SetGenerator(generator);
    }
  writer->SetMetaData(meta);
  writer->SetFilePrefix(directoryName);
  writer->SetFilePattern("%s/IM-0001-%04.4d.dcm");
  writer->TimeAsVectorOn();
  if (reader)
    {
    if (reader->GetTimeDimension() > 1)
      {
      writer->SetTimeDimension(reader->GetTimeDimension());
      writer->SetTimeSpacing(reader->GetTimeSpacing());
      }
    if (reader->GetRescaleSlope() > 0)
      {
      writer->SetRescaleSlope(reader->GetRescaleSlope());
      writer->SetRescaleIntercept(reader->GetRescaleIntercept());
      }
    writer->SetMemoryRowOrder(reader->GetMemoryRowOrder());
    }
  writer->SET_INPUT_DATA(data);
  writer->SetPatientMatrix(matrix);
  writer->Write();
}

#else

vtkDICOMImageReader *ReadDICOMImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *directoryName,
  int coordSystem)
{
  // read the image
  vtkDICOMImageReader *reader = vtkDICOMImageReader::New();

  reader->SetDirectoryName(directoryName);
  reader->Update();
  if (reader->GetErrorCode())
    {
    exit(1);
    }

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

  return reader;
}
#endif

vtkMINCImageReader *ReadMINCImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *fileName,
  int coordSystem)
{
  // read the image
  vtkMINCImageReader *reader = vtkMINCImageReader::New();

  reader->SetFileName(fileName);
  reader->Update();
  if (reader->GetErrorCode())
    {
    exit(1);
    }

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

  return reader;
}

void WriteMINCImage(
  vtkImageReader2 *vtkNotUsed(sourceReader),
  vtkImageReader2 *vtkNotUsed(targetReader),
  vtkImageData *data, vtkMatrix4x4 *vtkNotUsed(matrix), const char *fileName,
  int vtkNotUsed(coordSystem))
{
  fprintf(stderr, "Writing MINC images is not supported yet, "
          "the output file will have incorrect information\n");
  vtkSmartPointer<vtkMINCImageWriter> writer =
    vtkSmartPointer<vtkMINCImageWriter>::New();
  writer->SetFileName(fileName);
  writer->SET_INPUT_DATA(data);
  // the input matrix must be converted
  //writer->SetDirectionCosines(matrix);
  writer->Write();
}

#ifdef AIRS_USE_NIFTI
vtkNIFTIReader *ReadNIFTIImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *fileName,
  int coordSystem)
{
  // read the image
  vtkNIFTIReader *reader = vtkNIFTIReader::New();

  reader->SetFileName(fileName);
  reader->Update();
  if (reader->GetErrorCode())
    {
    exit(1);
    }

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
  if (reader->GetQFormMatrix())
    {
    vtkMatrix4x4::DeepCopy(nMatrix, reader->GetQFormMatrix());
    }
  else if (reader->GetSFormMatrix())
    {
    vtkMatrix4x4::DeepCopy(nMatrix, reader->GetSFormMatrix());
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

  return reader;
}

void WriteNIFTIImage(
  vtkImageReader2 *vtkNotUsed(sourceReader), vtkImageReader2 *targetReader,
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *fileName,
  int vtkNotUsed(coordSystem))
{
  vtkNIFTIReader *reader = vtkNIFTIReader::SafeDownCast(targetReader);

  vtkSmartPointer<vtkNIFTIWriter> writer =
    vtkSmartPointer<vtkNIFTIWriter>::New();
  if (reader)
    {
    writer->SetNIFTIHeader(reader->GetNIFTIHeader());
    if (reader->GetTimeDimension() > 1)
      {
      writer->SetTimeDimension(reader->GetTimeDimension());
      writer->SetTimeSpacing(reader->GetTimeSpacing());
      }
    if (reader->GetQFac() < 0)
      {
      writer->SetQFac(-1.0);
      }
    }
  writer->SET_INPUT_DATA(data);
  writer->SetQFormMatrix(matrix);
  writer->SetSFormMatrix(matrix);
  writer->SetFileName(fileName);
  writer->Write();
}

#endif /* AIRS_USE_NIFTI */

vtkImageReader2 *ReadImage(
  vtkImageData *image, vtkMatrix4x4 *matrix,
  const char *filename, int coordSystem)
{
  int t = GuessFileType(filename);

  if (t == MINCImage)
    {
    return ReadMINCImage(image, matrix, filename, coordSystem);
    }
  else if (t == NIFTIImage)
    {
#ifdef AIRS_USE_NIFTI
    return ReadNIFTIImage(image, matrix, filename, coordSystem);
#else
    fprintf(stderr, "NIFTI files are not supported.\n");
    exit(1);
#endif
    }

  return ReadDICOMImage(image, matrix, filename, coordSystem);
}

int CoordSystem(const char *filename)
{
  int t = GuessFileType(filename);

  if (t == MINCImage || t == NIFTIImage)
    {
    return NIFTICoords;
    }

  return DICOMCoords;
}

void WriteImage(
  vtkImageReader2 *sourceReader, vtkImageReader2 *targetReader,
  vtkImageData *image, vtkMatrix4x4 *matrix,
  const char *filename, int coordSystem)
{
  int t = GuessFileType(filename);

  if (t == MINCImage)
    {
    WriteMINCImage(
      sourceReader, targetReader, image, matrix, filename, coordSystem);
    }
  else if (t == NIFTIImage)
    {
#ifdef AIRS_USE_NIFTI
    WriteNIFTIImage(
      sourceReader, targetReader, image, matrix, filename, coordSystem);
#else
    fprintf(stderr, "NIFTI files are not supported.\n");
    exit(1);
#endif
    }
  else
    {
#ifdef AIRS_USE_DICOM
    WriteDICOMImage(
      sourceReader, targetReader, image, matrix, filename, coordSystem);
#else
    fprintf(stderr, "Writing DICOM files is not supported.\n");
    exit(1);
#endif
    }
}


void WriteMesh(vtkPolyData *mesh, vtkMatrix4x4 *matrix,
  const char *filename)
{
  int t = GuessFileType(filename);

  vtkSmartPointer<vtkTransform> transform =
    vtkSmartPointer<vtkTransform>::New();
  if (matrix)
    {
    transform->Concatenate(matrix);
    }

  vtkSmartPointer<vtkTransformPolyDataFilter> filter =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  filter->SET_INPUT_DATA(mesh);
  filter->SetTransform(transform);
  filter->Update();

  if (t == STLSurface)
    {
    vtkSmartPointer<vtkSTLWriter> writer =
      vtkSmartPointer<vtkSTLWriter>::New();
    writer->SET_INPUT_DATA(filter->GetOutput());
    writer->SetFileName(filename);
    writer->Write();
    }
  else if (t == OBJSurface)
    {
    vtkSmartPointer<vtkMNIObjectWriter> writer =
      vtkSmartPointer<vtkMNIObjectWriter>::New();
    writer->SET_INPUT_DATA(filter->GetOutput());
    writer->SetFileName(filename);
    writer->Write();
    }
  else if (t == VTKSurface)
    {
    vtkSmartPointer<vtkPolyDataWriter> writer =
      vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SET_INPUT_DATA(filter->GetOutput());
    writer->SetFileName(filename);
    writer->Write();
    }
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

// a class to look for errors when reading transforms.
class ErrorObserver : public vtkCommand
{
public:
  static ErrorObserver *New() { return new ErrorObserver; }
  vtkTypeMacro(ErrorObserver, vtkCommand);
  virtual void Execute(vtkObject *o, unsigned long eventId, void *callData);
};

void ErrorObserver::Execute(
  vtkObject *, unsigned long, void *callData)
{
  if (callData)
    {
    fprintf(stderr, "%s\n", static_cast<char *>(callData));
    }
  exit(1);
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

void ComputeRange(vtkImageData *image, double range[2])
{
  // compute the range within a cylinder that is slightly smaller than
  // the image bounds (the idea is to capture only the reconstructed
  // portion of a CT image).
  double spacing[3];
  double origin[3];
  int extent[6];
  double bounds[6];
  image->GetSpacing(spacing);
  image->GetOrigin(origin);
  image->GetExtent(extent);

  for (int i = 0; i < 3; ++i)
    {
    double b1 = extent[2*i]*spacing[i] + origin[i];
    double b2 = extent[2*i+1]*spacing[i] + origin[i];
    bounds[2*i] = (b1 < b2 ? b1 : b2);
    bounds[2*i+1] = (b1 < b2 ? b2 : b1);
    spacing[i] = fabs(spacing[i]);
    origin[i] = bounds[2*i];
    // reduce bounds by 2% in X and Y for use in cylinder generation
    double bl = (i == 2 ? 0.0 : 0.01*(bounds[2*i+1] - bounds[2*i]));
    bounds[2*i] += bl;
    bounds[2*i+1] -= bl;
    }

  // extract just the reconstructed portion of CT image
  vtkSmartPointer<vtkROIStencilSource> cylinder =
    vtkSmartPointer<vtkROIStencilSource>::New();

  cylinder->SetShapeToCylinderZ();
  cylinder->SetInformationInput(image);
  cylinder->SetBounds(bounds);
  cylinder->Update();

  // get the range within the cylinder
  vtkSmartPointer<vtkImageHistogramStatistics> rangeFinder =
    vtkSmartPointer<vtkImageHistogramStatistics>::New();

  rangeFinder->GetAutoRangePercentiles(range);
  rangeFinder->SET_INPUT_DATA(image);
  rangeFinder->SET_STENCIL_DATA(cylinder->GetOutput());
  rangeFinder->Update();

  rangeFinder->GetAutoRange(range);
}

};

struct skullstrip_options
{
  double bt;           // -t --threshold
  double d1;           // --d1
  double d2;           // --d2
  double t;            // --tesselations
  double rmin;         // --rmin
  double rmax;         // --rmax
  int n;               // -N --iterations
  int display;         // -d --display
  int silent;          // -s --silent
  int coords;          // -C --coords
  const char *surface; // -o (output mesh)
  const char *output;  // -o (output image)
  const char *screenshot; // -j (output screenshot)
  const char *source;
};

void skullstrip_initialize_options(skullstrip_options *options)
{
  options->bt = 0.5;
  options->d1 = 20.0;
  options->d2 = 10.0;
  options->rmin = 3.0;
  options->rmax = 10.0;
  options->n = 1000;
  options->t = 4;
  options->display = 0;
  options->silent = 0;
  options->coords = NativeCoords;
  options->screenshot = NULL;
  options->output = NULL;
  options->surface = NULL;
  options->source = NULL;
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

  bool match = false;
  for (const char **t = possib; *t != 0; t++)
    {
    if (strcmp(*t, arg) == 0)
      {
      match = true;
      break;
      }
    }

  if (!match)
    {
    fprintf(stderr, "Incorrect value for option \"%s\": %s\n",
            op, arg);
    fprintf(stderr, "Allowed values:");
    for (const char **u = possib; *u != 0; u++)
      {
      fprintf(stderr, "%s", *u);
      }
    fprintf(stderr, "\n");
    exit(1);
    }

  return arg;
}

void skullstrip_show_usage(FILE *fp, const char *command)
{
  const char *cp = command + strlen(command);
  while (cp > command && cp[-1] != '/' && cp[-1] != '\\') { --cp; }

  fprintf(fp,
    "Usage: %s [options] -o <output> <source image>\n", cp);
  fprintf(fp, "\n");
  fprintf(fp,
    "For more information, type \"%s --help\"\n\n", command);
}

void skullstrip_show_help(FILE *fp, const char *command)
{
  const char *cp = command + strlen(command);
  while (cp > command && cp[-1] != '/' && cp[-1] != '\\') { --cp; }

  fprintf(fp,
    "Usage: %s [options] -o <output> <source image>\n", cp);
  fprintf(fp,
    "\n"
    "Written by David Gobbi <dgobbi@ucalgary.ca> at CIPAC.  Version 0.1.\n"
    "\n"
    "This program takes a 3D structural head MRI as input, and produces an\n"
    "output image that consists of only the brain, with everything colored\n"
    "black.  It uses the BET (Brain Extraction Technique) algorithm, and\n"
    "has all of the strengths and weaknesses of that algorithm: it is very\n"
    "fast, but undersegments the brain in certain regions, in particular\n"
    "around the orbits.  The parameters of the algorithm can be adjusted to\n"
    "improve the results (within limits).\n"
    "\n"
    "The \"-o\" option allows you to specify either an output image, or an\n"
    "output mesh (geodesic surface).  For an output image, supply a file\n"
    "with a suffix of .nii (NIFTI), .mnc (MINC), or supply the name of an\n"
    "existing directory (for DICOM).  For an output mesh, provide a file\n"
    "with a suffix of .stl, .vtk, or .obj.\n"
    "\n");
  fprintf(fp,
    " --threshold (default: 0.5)\n"
    "\n"
    "    The threshold that will be used to detect the edge of the brain,\n"
    "    as a fraction of the 98th percentile of pixel values.  The 98th\n"
    "    percentile is an estimate of the average of the brightest\n"
    "    pixels in the image, excluding outliers.\n"
    "\n"
    "  --d1 (default: 20.0)\n"
    "\n"
    "    Distance (in mm) that is searched for minimum intensity values\n"
    "    during segmentation of the brain surface.\n"
    "\n"
    "  --d2 (default: 10.0)\n"
    "\n"
    "    Distance (in mm) that is searched for maximum intensity values\n"
    "    during segmentation of the brain surface.\n"
    "\n"
    "  --rmin (default: 3.0)\n"
    "\n"
    "    Minimum radius of curvature.\n"
    "\n"
    "  --rmin (default: 10.0)\n"
    "\n"
    "    Maxiumum radius of curvature at which smoothing will be applied.\n"
    "\n"
    " -N --iterations (default: 500x500x500)\n"
    "\n"
    "    Set the maximum number of iterations per stage.  Set this to zero\n"
    "    if you want to use the initial transform as-is.\n"
    "\n"
    " -C --coords (default: guess from file type)\n"
    "                 DICOM     LPS\n"
    "                 NIFTI     RAS\n"
    "                 MINC      RAS\n"
    "\n"
    "    The DICOM standard defines a patient coordinate system where the\n"
    "    x-axis points left, the y-axis points towards the back, and the\n"
    "    z-axis points towards the head.  NIFTI and MINC use a coordinate\n"
    "    system where x points right and y points towards the front.  The\n"
    "    mesh that is written by the \"-o\" option will be in the chosen\n"
    "    coordinate system.\n"
    "\n"
    " -d --display      (default: off)\n"
    "\n"
    "    Display the result.\n"
    "\n"
    " -s --silent       (default: off)\n"
    "\n"
    "    Do not print information to the console during the stripping.\n"
    "    This is useful when running in batch mode.  Error messages will\n"
    "    still be printed.\n"
    "\n"
    " -j --screenshot <file>\n"
    "\n"
    "    Write a screenshot as a png, jpeg, or tiff file.  This is useful\n"
    "    when performing stripping in a batch file in order to provide\n"
    "    a simple means of visually assessing the results retrospectively.\n"
    "\n"
    " -o <file>\n"
    "\n"
    "    Provide an output file, either an image or a mesh.  If both an\n"
    "    image and a mesh are to be written as output, then use the -o\n"
    "    option twice, once for each.\n"
    "\n"
    "Reference for algorithm:\n"
    "\n"
    "    Stephen M. Smith, \"Fast Robust Automated Brain Extraction,\" Human\n"
    "    Brain Mapping 17:143-155, 2002\n"
    "\n");
}

int skullstrip_read_options(
  int argc, char *argv[], skullstrip_options *options)
{
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
      int t = GuessFileType(arg);

      if (t <= LastImageType)
        {
        if (options->source == 0)
          {
          options->source = arg;
          }
        else
          {
          fprintf(stderr, "Too many input images listed on command line\n");
          exit(1);
          }
        }
      }
    else
      {
      if (strcmp(arg, "-h") == 0 ||
          strcmp(arg, "--help") == 0)
        {
        skullstrip_show_help(stdout, argv[0]);
        exit(0);
        }
      else if (strcmp(arg, "--threshold") == 0)
        {
        arg = check_next_arg(argc, argv, &argi, 0);
        options->bt = strtod(arg, const_cast<char **>(&arg));
        }
      else if (strcmp(arg, "--d1") == 0)
        {
        arg = check_next_arg(argc, argv, &argi, 0);
        options->d1 = strtod(arg, const_cast<char **>(&arg));
        }
      else if (strcmp(arg, "--d2") == 0)
        {
        arg = check_next_arg(argc, argv, &argi, 0);
        options->d2 = strtod(arg, const_cast<char **>(&arg));
        }
      else if (strcmp(arg, "--rmin") == 0)
        {
        arg = check_next_arg(argc, argv, &argi, 0);
        options->rmin = strtod(arg, const_cast<char **>(&arg));
        }
      else if (strcmp(arg, "--rmax") == 0)
        {
        arg = check_next_arg(argc, argv, &argi, 0);
        options->rmax = strtod(arg, const_cast<char **>(&arg));
        }
      else if (strcmp(arg, "--iterations") == 0 ||
               strcmp(arg, "-N") == 0)
        {
        arg = check_next_arg(argc, argv, &argi, 0);
        options->n = strtoul(arg, const_cast<char **>(&arg), 0);
        }
      else if (strcmp(arg, "--tesselations") == 0)
        {
        arg = check_next_arg(argc, argv, &argi, 0);
        options->t = strtoul(arg, const_cast<char **>(&arg), 0);
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
      else if (strcmp(arg, "-d") == 0 ||
               strcmp(arg, "--display") == 0)
        {
        options->display = 1;
        }
      else if (strcmp(arg, "-s") == 0 ||
               strcmp(arg, "--silent") == 0)
        {
        options->silent = 1;
        }
      else if (strcmp(arg, "-j") == 0 ||
               strcmp(arg, "--screenshot") == 0)
        {
        arg = check_next_arg(argc, argv, &argi, 0);
        options->screenshot = arg;
        }
      else if (strcmp(arg, "-o") == 0)
        {
        arg = check_next_arg(argc, argv, &argi, 0);
        int t = GuessFileType(arg);
        if (t <= LastImageType)
          {
          if (options->output)
            {
            fprintf(stderr, "Too many -o options specified!\n");
            skullstrip_show_usage(stderr, argv[0]);
            }
          options->output = arg;
          }
        else if (t <= LastSurfaceType)
          {
          if (options->surface)
            {
            fprintf(stderr, "Too many -o options specified!\n");
            skullstrip_show_usage(stderr, argv[0]);
            }
          options->surface = arg;
          }
        }
      else
        {
        fprintf(stderr, "Unrecognized option \"%s\"\n", arg);
        skullstrip_show_usage(stderr, argv[0]);
        exit(1);
        }
      }
    }

  return 1;
}

int main(int argc, char *argv[])
{
  skullstrip_options options;
  skullstrip_initialize_options(&options);
  skullstrip_read_options(argc, argv, &options);

  // -------------------------------------------------------
  // the files
  const char *imagefile = options.output;
  const char *meshfile = options.surface;
  const char *sourcefile = options.source;
  bool display = (options.display != 0 ||
                  options.screenshot != 0);

  if (!sourcefile)
    {
    skullstrip_show_usage(stderr, argv[0]);
    return 1;
    }

  // -------------------------------------------------------
  // load the images

  if (options.coords == NativeCoords)
    {
    int ic = CoordSystem(sourcefile);

    if (ic == DICOMCoords)
      {
      options.coords = DICOMCoords;
      }
    else
      {
      options.coords = NIFTICoords;
      }
    }

  if (!options.silent)
    {
    cout << "Reading source image: " << sourcefile << endl;
    }

  vtkSmartPointer<vtkImageData> sourceImage =
    vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkMatrix4x4> sourceMatrix =
    vtkSmartPointer<vtkMatrix4x4>::New();
  vtkSmartPointer<vtkImageReader2> sourceReader =
    ReadImage(sourceImage, sourceMatrix, sourcefile, options.coords);
  sourceReader->Delete();

  if (!options.silent)
    {
    if (options.coords == DICOMCoords)
      {
      cout << "Using DICOM patient coords." << endl;;
      }
    else
      {
      cout << "Using NIFTI (or MINC) world coords." << endl;
      }
    }

  // -------------------------------------------------------
  // make a timer
  vtkSmartPointer<vtkTimerLog> timer =
    vtkSmartPointer<vtkTimerLog>::New();
  double startTime = timer->GetUniversalTime();

  // -------------------------------------------------------
  // find the direction that is "superior"
  double maxval = 0;
  int maxidx = 0;
  for (int idx = 0; idx < 3; idx++)
    {
    double vec[4] = { 0.0, 0.0, 0.0, 0.0 };
    vec[idx] = 1.0;
    sourceMatrix->MultiplyPoint(vec, vec);
    if (fabs(vec[2]) > fabs(maxval))
      {
      maxval = vec[2];
      maxidx = idx;
      }
    }

  // set the Z extent to a set fraction of the other extents
  int brainExtent[6];
  sourceImage->GetExtent(brainExtent);
  int idx1 = (maxidx + 1) % 3;
  int idx2 = (idx1 + 1) % 3;
  int size1 = brainExtent[2*idx1 + 1] - brainExtent[2*idx1] + 1;
  int size2 = brainExtent[2*idx2 + 1] - brainExtent[2*idx2] + 1;
  int maxsize = brainExtent[2*maxidx + 1] - brainExtent[2*maxidx] + 1;
  int size3 = vtkMath::Floor(sqrt(0.7*0.7*size1*size2));
  // make sure computed size is not more than full image size
  if (size3 > maxsize)
    {
    size3 = maxsize;
    }
  if (maxval < 0)
    {
    brainExtent[2*maxidx + 1] -= (maxsize - size3);
    }
  else
    {
    brainExtent[2*maxidx] += (maxsize - size3);
    }

  // -------------------------------------------------------
  // strip the image

  vtkSmartPointer<vtkImageMRIBrainExtractor> stripper =
    vtkSmartPointer<vtkImageMRIBrainExtractor>::New();
  stripper->SET_INPUT_DATA(sourceImage);
  stripper->SetRMin(options.rmin);
  stripper->SetRMax(options.rmax);
  stripper->SetD1(options.d1);
  stripper->SetD2(options.d2);
  stripper->SetBT(options.bt);
  stripper->SetNumberOfIterations(options.n);
  stripper->SetNumberOfTessellations(options.t);
  stripper->SetBrainExtent(brainExtent);
  stripper->Update();

  double lastTime = timer->GetUniversalTime();

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

  sourceMapper->SET_INPUT_DATA(sourceImage);
  sourceMapper->SliceAtFocalPointOn();
  sourceMapper->SliceFacesCameraOn();
  sourceMapper->ResampleToScreenPixelsOff();

  double sourceRange[2];
  sourceImage->GetScalarRange(sourceRange);
  ComputeRange(sourceImage, sourceRange);

  sourceProperty->SetInterpolationTypeToLinear();
  sourceProperty->SetColorWindow((sourceRange[1]-sourceRange[0]));
  sourceProperty->SetColorLevel(0.5*(sourceRange[0]+sourceRange[1]));

  sourceActor->SetMapper(sourceMapper);
  sourceActor->SetProperty(sourceProperty);
  sourceActor->SetUserMatrix(sourceMatrix);

  vtkSmartPointer<vtkImageSlice> brainActor =
    vtkSmartPointer<vtkImageSlice>::New();
  vtkSmartPointer<vtkImageResliceMapper> brainMapper =
    vtkSmartPointer<vtkImageResliceMapper>::New();
  vtkSmartPointer<vtkImageProperty> brainProperty =
    vtkSmartPointer<vtkImageProperty>::New();
  vtkSmartPointer<vtkLookupTable> brainTable =
    vtkSmartPointer<vtkLookupTable>::New();

  brainMapper->SET_INPUT_DATA(stripper->GetOutput());
  brainMapper->SliceAtFocalPointOn();
  brainMapper->SliceFacesCameraOn();
  brainMapper->ResampleToScreenPixelsOff();

  brainTable->SetRampToLinear();
  brainTable->SetHueRange(0.0, 0.0);
  brainTable->SetSaturationRange(1.0, 1.0);
  brainTable->SetValueRange(0.0, 1.0);
  brainTable->Build();
  brainTable->SetTableValue(0, 0.0, 0.0, 0.0, 0.0);

  brainProperty->SetInterpolationTypeToLinear();
  brainProperty->SetColorWindow((sourceRange[1]-sourceRange[0]));
  brainProperty->SetColorLevel(0.5*(sourceRange[0]+sourceRange[1]));
  brainProperty->SetLookupTable(brainTable);

  brainActor->SetMapper(brainMapper);
  brainActor->SetProperty(brainProperty);
  brainActor->SetUserMatrix(sourceMatrix);

  vtkSmartPointer<vtkImageStack> imageStack =
    vtkSmartPointer<vtkImageStack>::New();
  imageStack->AddImage(sourceActor);
  imageStack->AddImage(brainActor);

  renderer->AddViewProp(imageStack);
  renderer->SetBackground(0,0,0);

  renderWindow->SetSize(512,512);

  double bounds[6], center[4], tspacing[3];
  int extent[6];
  sourceImage->GetBounds(bounds);
  sourceImage->GetExtent(extent);
  sourceImage->GetSpacing(tspacing);
  center[0] = 0.5*(bounds[0] + bounds[1]);
  center[1] = 0.5*(bounds[2] + bounds[3]);
  center[2] = 0.5*(bounds[4] + bounds[5]);
  center[3] = 1.0;
  sourceMatrix->MultiplyPoint(center, center);

  vtkCamera *camera = renderer->GetActiveCamera();
  renderer->ResetCamera();
  camera->SetFocalPoint(center);
  camera->ParallelProjectionOn();
  camera->SetParallelScale(0.5*(bounds[3] - bounds[2]));
  SetViewFromMatrix(renderer, istyle, sourceMatrix, options.coords);
  renderer->ResetCameraClippingRange();

  if (display)
    {
    renderWindow->Render();
    }

  if (!options.silent)
    {
    cout << "stripping took " << (lastTime - startTime) << "s" << endl;
    }

  // -------------------------------------------------------
  // capture a screen shot
  if (options.screenshot)
    {
    WriteScreenshot(renderWindow, options.screenshot);
    }

  // -------------------------------------------------------
  // write the output file
  if (imagefile)
    {
    if (!options.silent)
      {
      cout << "Writing stripped image: " << imagefile << endl;
      }

    WriteImage(sourceReader, sourceReader,
      stripper->GetOutput(), sourceMatrix, imagefile, options.coords);
    }
  if (meshfile)
    {
    if (!options.silent)
      {
      cout << "Writing brain surface: " << meshfile << endl;
      }

    WriteMesh(stripper->GetBrainMesh(), sourceMatrix, meshfile);
    }

  if (!options.silent)
    {
    cout << "Done!" << endl;
    }

  // -------------------------------------------------------
  // allow user to interact

  if (options.display)
    {
    interactor->Start();
    }

  return 0;
}
