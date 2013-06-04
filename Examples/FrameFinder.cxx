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


#include <vtkSmartPointer.h>

#include <vtkImageReslice.h>
#include <vtkImageResize.h>
#include <vtkImageSincInterpolator.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkTubeFilter.h>

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

#include <vtkFrameFinder.h>

// internal methods for reading images, these methods read the image
// into the specified data object and also provide a matrix for converting
// the data coordinates into patient coordinates.
namespace {

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

  // -------------------------------------------------------
  // make the frame
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(16);

  vtkSmartPointer<vtkCellArray> cells =
    vtkSmartPointer<vtkCellArray>::New();

  static double leksellPoints[4][4][3] = {
    { { 196.0, 40.0, 160.0 }, { 196.0, 40.0, 40.0 },
      { 196.0, 160.0, 160.0 }, {196.0, 160.0, 40.0 } },
    { { 4.0, 40.0, 160.0 }, { 4.0, 40.0, 40.0 },
      { 4.0, 160.0, 160.0 }, { 4.0, 160.0, 40.0 } },
    { { 40.0, 217.5, 160.0 }, { 40.0, 217.5, 40.0 },
      { 160.0, 217.5, 160.0 }, { 160.0, 217.5, 40.0 } },
    { { 40.0, -17.5, 160.0 }, { 40.0, -17.5, 40.0 },
      { 160.0, -17.5, 160.0 }, { 160.0, -17.5, 40.0 } } };

  for (int i = 0; i < 4; i++)
    {
    double (*p)[3] = leksellPoints[i];

    cells->InsertNextCell(4);
    for (int j = 0; j < 4; j++)
      {
      int ptIdx = i*4 + j;
      points->SetPoint(ptIdx, p[j]);
      cells->InsertCellPoint(ptIdx);
      }
    }

  vtkSmartPointer<vtkPolyData> frameData =
    vtkSmartPointer<vtkPolyData>::New();
  frameData->SetPoints(points);
  frameData->SetLines(cells);

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

  sourceMapper->SetInput(sourceImage);
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

  vtkSmartPointer<vtkTubeFilter> tube =
    vtkSmartPointer<vtkTubeFilter>::New();
  tube->SetInput(frameData);
  tube->SetRadius(0.5);
  tube->SetNumberOfSides(10);

  vtkSmartPointer<vtkDataSetMapper> targetMapper =
    vtkSmartPointer<vtkDataSetMapper>::New();
  //targetMapper->SetInputConnection(tube->GetOutputPort());
  targetMapper->SetInput(frameData);

  vtkSmartPointer<vtkActor> targetActor =
    vtkSmartPointer<vtkActor>::New();
  targetActor->SetMapper(targetMapper);
  targetActor->SetUserMatrix(targetMatrix);
  targetActor->GetProperty()->SetAmbient(0.6);
  targetActor->GetProperty()->SetColor(1.0, 1.0, 0.0);
  targetActor->GetProperty()->SetOpacity(0.5);

  renderer->AddViewProp(sourceActor);
  renderer->AddViewProp(targetActor);
  renderer->SetBackground(0.2,0.2,0.2);

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
  frameFinder->SetInput(sourceImage);

  // -------------------------------------------------------
  // make a timer
  vtkSmartPointer<vtkTimerLog> timer =
    vtkSmartPointer<vtkTimerLog>::New();
  double startTime = timer->GetUniversalTime();
  double lastTime = startTime;

  // -------------------------------------------------------
  // do the frame finding

  frameFinder->Update();

  vtkSmartPointer<vtkMatrix4x4> frameMatrix =
    vtkSmartPointer<vtkMatrix4x4>::New();
  vtkMatrix4x4::Multiply4x4(
    leksellToDICOM, frameFinder->GetImageToFrameMatrix(), frameMatrix);
  sourceActor->SetUserMatrix(frameMatrix);

  vtkSmartPointer<vtkDataSetMapper> frameMapper =
    vtkSmartPointer<vtkDataSetMapper>::New();
  frameMapper->SetInputConnection(frameFinder->GetOutputPort());

  vtkSmartPointer<vtkActor> frameActor =
    vtkSmartPointer<vtkActor>::New();
  frameActor->SetUserMatrix(frameMatrix);
  frameActor->SetMapper(frameMapper);
  frameActor->GetProperty()->SetAmbient(0.6);
  frameActor->GetProperty()->SetColor(1.0, 0.0, 1.0);

  renderer->AddViewProp(frameActor);
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

  // -------------------------------------------------------
  // allow user to interact

  interactor->Start();

  return 1;
}
