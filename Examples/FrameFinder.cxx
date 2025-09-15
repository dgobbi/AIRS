/*=========================================================================

Program:   Atamai Image Registration and Segmentation
Module:    FrameFinder.cxx

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

#include "AIRSConfig.h"

#include <vtkSmartPointer.h>

#include <vtkImageReslice.h>
#include <vtkImageResize.h>
#include <vtkImageSincInterpolator.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkMath.h>
#include <vtkLookupTable.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkTubeFilter.h>
#include <vtkSphereSource.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkGlyph3D.h>
#include <vtkStringArray.h>

#include <vtkMINCImageReader.h>
#include <vtkDICOMImageReader.h>
#include <vtkMNITransformWriter.h>

#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageSlice.h>
#include <vtkImageStack.h>
#include <vtkImageResliceMapper.h>
#include <vtkImageProperty.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

#include <vtkTimerLog.h>
#include <vtkVersion.h>

#include <vtkFrameFinder.h>

// optional readers
#ifdef AIRS_USE_DICOM
#include <vtkDICOMReader.h>
#include <vtkDICOMSorter.h>
#include <vtkDICOMMetaData.h>
#include <vtkGlobFileNames.h>
#endif

#include <vtksys/SystemTools.hxx>
#include <string>

// A macro to assist VTK 5 backwards compatibility
#if VTK_MAJOR_VERSION >= 6
#define SET_INPUT_DATA SetInputData
#else
#define SET_INPUT_DATA SetInput
#endif

// internal methods for reading images, these methods read the image
// into the specified data object and also provide a matrix for converting
// the data coordinates into patient coordinates.
namespace {

#ifdef AIRS_USE_DICOM
void ReadDICOMImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *directoryName)
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
  vtkSmartPointer<vtkDICOMReader> reader =
    vtkSmartPointer<vtkDICOMReader>::New();
  reader->SetFileNames(sorter->GetFileNamesForSeries(0));
  reader->SetMemoryRowOrderToFileNative();

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
}

#else /* AIRS_USE_DICOM */

void ReadDICOMImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *directoryName)
{
  // read the image
  vtkSmartPointer<vtkDICOMImageReader> reader =
    vtkSmartPointer<vtkDICOMImageReader>::New();

  reader->SetDirectoryName(directoryName);
  reader->Update();

  // the reader flips the image and reverses the ordering, so undo these
  vtkSmartPointer<vtkImageReslice> flip =
    vtkSmartPointer<vtkImageReslice>::New();

  flip->SetInputConnection(reader->GetOutputPort());
  flip->SetResliceAxesDirectionCosines(
    1,0,0, 0,-1,0, 0,0,-1);
  flip->Update();

  vtkImageData *image = flip->GetOutput();

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
  matrix->Modified();
}

#endif /* AIRS_USE_DICOM */

void ReadMINCImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *fileName)
{
  // read the image
  vtkSmartPointer<vtkMINCImageReader> reader =
    vtkSmartPointer<vtkMINCImageReader>::New();

  reader->SetFileName(fileName);
  reader->Update();

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

  vtkImageData *image = flip->GetOutput();

  // get the data
  data->CopyStructure(image);
  data->GetPointData()->PassData(image->GetPointData());

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

void SetViewFromMatrix(
  vtkRenderer *renderer,
  vtkInteractorStyleImage *istyle,
  vtkMatrix4x4 *matrix)
{
  istyle->SetCurrentRenderer(renderer);

  // This view assumes the data uses the DICOM Patient Coordinate System.
  // It provides a right-is-left view of axial and coronal images
  double viewRight[4] = { 1.0, 0.0, 0.0, 0.0 };
  double viewUp[4] = { 0.0, -1.0, 0.0, 0.0 };

  matrix->MultiplyPoint(viewRight, viewRight);
  matrix->MultiplyPoint(viewUp, viewUp);

  istyle->SetImageOrientation(viewRight, viewUp);
}

};

void printUsage(const char *cmdname)
{
    cout << "Usage 1: " << cmdname << " --nodisplay -o output.xfm file.mnc"
         << endl;
    cout << "Usage 2: " << cmdname << " --nodisplay -o output.xfm dicomdir/"
         << endl;
}

int main (int argc, char *argv[])
{
  if (argc < 2)
  {
    printUsage(argv[0]);
    return EXIT_FAILURE;
  }

  // -------------------------------------------------------
  // the files
  int argi = 1;
  const char *xfmfile = NULL;
  const char *sourcefile;
  bool display = true;

  if (strcmp(argv[argi], "--nodisplay") == 0)
  {
    display = false;
    argi++;
  }
  if (strcmp(argv[argi], "-o") == 0)
  {
    if (argc <= argi + 1)
    {
      cerr << argv[0] << " : missing .xfm file after -o\n" << endl;
      return EXIT_FAILURE;
    }
    xfmfile = argv[argi + 1];
    argi += 2;
    size_t m = strlen(xfmfile);
    if (m < 4 || strcmp(&xfmfile[m-4], ".xfm") != 0)
    {
      cerr << argv[0] << " : transform file must end in .xfm\n" << endl;
      return EXIT_FAILURE;
    }
  }

  if (argc <= argi)
  {
    printUsage(argv[0]);
    return EXIT_FAILURE;
  }
  sourcefile = argv[argi];

  // -------------------------------------------------------
  // load the images

  int n = 0;

  // Read the source image
  vtkSmartPointer<vtkImageData> sourceImage =
    vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkMatrix4x4> sourceMatrix =
    vtkSmartPointer<vtkMatrix4x4>::New();
  n = strlen(sourcefile);
  if (n > 4 && strcmp(&sourcefile[n-4], ".mnc") == 0)
  {
    ReadMINCImage(sourceImage, sourceMatrix, sourcefile);
  }
  else
  {
    ReadDICOMImage(sourceImage, sourceMatrix, sourcefile);
  }

  static double leksellToDICOM16[16] = {
     1.0,  0.0,  0.0, -100.0,
     0.0, -1.0,  0.0,  100.0,
     0.0,  0.0, -1.0,  100.0,
     0.0,  0.0,  0.0,    1.0 };

  vtkSmartPointer<vtkMatrix4x4> leksellToDICOM =
    vtkMatrix4x4::New();
  leksellToDICOM->DeepCopy(leksellToDICOM16);

  vtkSmartPointer<vtkMatrix4x4> targetMatrix =
    vtkMatrix4x4::New();
  targetMatrix->DeepCopy(leksellToDICOM);

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

  istyle->SetInteractionModeToImage3D(); //Slicing();
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
  //sourceMapper->SliceFacesCameraOn();
  sourceMapper->ResampleToScreenPixelsOff();

  double sourceRange[2];
  sourceImage->GetScalarRange(sourceRange);
  sourceProperty->SetInterpolationTypeToLinear();
  sourceProperty->SetColorWindow((sourceRange[1]-sourceRange[0]));
  sourceProperty->SetColorLevel(0.5*(sourceRange[0]+sourceRange[1]));

  sourceActor->SetMapper(sourceMapper);
  sourceActor->SetProperty(sourceProperty);
  sourceActor->SetUserMatrix(sourceMatrix);

  renderer->AddViewProp(sourceActor);
  renderer->SetBackground(0.2,0.2,0.2);
  //renderer->SetBackground(0.0,0.0,0.0);

  renderWindow->SetSize(720,720);

  //double bounds[6] = {-10.0, 210.0, -10.0, 210.0, -10.0, 210.0 };
  double center[4] = {100.0, 100.0, 100.0, 1.0};
  targetMatrix->MultiplyPoint(center, center);

  vtkCamera *camera = renderer->GetActiveCamera();
  renderer->ResetCamera();
  camera->SetFocalPoint(center);
  camera->ParallelProjectionOn();
  camera->SetParallelScale(132);
  SetViewFromMatrix(renderer, istyle, sourceMatrix);
  renderer->ResetCameraClippingRange();

  if (display)
  {
    renderWindow->Render();
  }

  // -------------------------------------------------------
  // prepare for frame finding

  vtkSmartPointer<vtkFrameFinder> frameFinder =
    vtkSmartPointer<vtkFrameFinder>::New();
  frameFinder->SetDICOMPatientMatrix(sourceMatrix);
  frameFinder->SET_INPUT_DATA(sourceImage);

  // -------------------------------------------------------
  // make a timer
  vtkSmartPointer<vtkTimerLog> timer =
    vtkSmartPointer<vtkTimerLog>::New();
  double startTime = timer->GetUniversalTime();
  double lastTime = startTime;

  // -------------------------------------------------------
  // do the frame finding

  frameFinder->Update();

  lastTime = timer->GetUniversalTime();

  vtkSmartPointer<vtkSphereSource> glyphSource =
    vtkSmartPointer<vtkSphereSource>::New();
  glyphSource->SetThetaResolution(21);
  //glyphSource->SetResolution(21);
  glyphSource->Update();

  vtkSmartPointer<vtkTransform> gtrans = vtkSmartPointer<vtkTransform>::New();
  //gtrans->RotateWXYZ(90, 1.0, 0.0, 0.0);

  vtkSmartPointer<vtkTransformPolyDataFilter> glyphRotate =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  glyphRotate->SetInputConnection(glyphSource->GetOutputPort());
  glyphRotate->SetTransform(gtrans);

  vtkSmartPointer<vtkGlyph3D> glypher =
    vtkSmartPointer<vtkGlyph3D>::New();
  glypher->SetSourceConnection(glyphRotate->GetOutputPort());
  glypher->SetInputConnection(frameFinder->GetOutputPort());
  glypher->SetScaleModeToScaleByVectorComponents();
  glypher->SetColorModeToColorByScalar();
  glypher->ScalingOn();
  glypher->OrientOff();
  glypher->ClampingOff();

  vtkSmartPointer<vtkMatrix4x4> frameMatrix =
    vtkSmartPointer<vtkMatrix4x4>::New();
  vtkMatrix4x4::Multiply4x4(
    leksellToDICOM, frameFinder->GetImageToFrameMatrix(), frameMatrix);
  sourceActor->SetUserMatrix(frameMatrix);

  vtkSmartPointer<vtkLookupTable> mtable =
    vtkSmartPointer<vtkLookupTable>::New();
  mtable->SetRampToLinear();
  mtable->SetSaturationRange(0.0, 0.0);
  mtable->SetValueRange(0.0, 1.0);
  mtable->Build();
  mtable->SetRange(sourceRange[0], sourceRange[1]*0.2);

  vtkSmartPointer<vtkDataSetMapper> frameMapper =
    vtkSmartPointer<vtkDataSetMapper>::New();
  frameMapper->SetInputConnection(glypher->GetOutputPort());
  frameMapper->SetColorModeToMapScalars();
  frameMapper->SetLookupTable(mtable);
  frameMapper->UseLookupTableScalarRangeOn();

  vtkSmartPointer<vtkActor> frameActor =
    vtkSmartPointer<vtkActor>::New();
  frameActor->SetUserMatrix(frameMatrix);
  frameActor->SetMapper(frameMapper);
  frameActor->GetProperty()->SetAmbient(0.1);
  frameActor->GetProperty()->SetDiffuse(1.0);
  //frameActor->GetProperty()->SetColor(1.0, 0.0, 1.0);

  renderer->AddViewProp(frameActor);

  vtkSmartPointer<vtkDataSetMapper> targetMapper =
    vtkSmartPointer<vtkDataSetMapper>::New();
  targetMapper->SetInputConnection(frameFinder->GetOutputPort(1));

  vtkSmartPointer<vtkActor> targetActor =
    vtkSmartPointer<vtkActor>::New();
  targetActor->SetMapper(targetMapper);
  targetActor->SetUserMatrix(targetMatrix);
  targetActor->GetProperty()->SetAmbient(0.6);
  targetActor->GetProperty()->SetColor(1.0, 1.0, 0.0);
  targetActor->GetProperty()->SetOpacity(0.5);

  renderer->AddViewProp(targetActor);

  if (display)
  {
    renderer->Render();
  }

/*
  // -------------------------------------------------------
  // write the output matrix
  if (xfmfile)
    {
    vtkSmartPointer<vtkMNITransformWriter> writer =
      vtkSmartPointer<vtkMNITransformWriter>::New();
    writer->SetFileName(xfmfile);
    writer->SetTransform(registration->GetTransform());
    writer->Update();
    }
*/

  cout << "Frame finding took " << (lastTime - startTime) << " seconds\n";

  // -------------------------------------------------------
  // allow user to interact

  interactor->Start();

  return 1;
}
