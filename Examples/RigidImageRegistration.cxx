/*=========================================================================

Program:   Atamai Image Registration and Segmentation
Module:    RigidImageRegistration.cxx

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

#include <vtkTimerLog.h>

#include <vtkImageRegistration.h>

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
    cout << "Usage 1: " << cmdname << " --nodisplay -o output.xfm source.mnc target.mnc"
         << endl;
    cout << "Usage 2: " << cmdname << " --nodisplay -o output.xfm dicomdir1/ dicomdir2/"
         << endl;
}

int main (int argc, char *argv[])
{
  if (argc < 3)
    {
    printUsage(argv[0]);
    return EXIT_FAILURE;
    }

  // -------------------------------------------------------
  // the files
  int argi = 1;
  const char *xfmfile = NULL;
  const char *sourcefile;
  const char *targetfile;
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

  if (argc <= argi + 1)
    {
    printUsage(argv[0]);
    return EXIT_FAILURE;
    }
  sourcefile = argv[argi];
  targetfile = argv[argi + 1];

  // -------------------------------------------------------
  // parameters for registration

  int interpolatorType = vtkImageRegistration::Linear;
  double transformTolerance = 0.1; // tolerance on transformation result
  int numberOfBins = 64; // for Mattes' mutual information
  double initialBlurFactor = 4.0;

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

  // Read the target image
  vtkSmartPointer<vtkImageData> targetImage =
    vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkMatrix4x4> targetMatrix =
    vtkSmartPointer<vtkMatrix4x4>::New();
  n = strlen(targetfile);
  if (n > 4 && strcmp(&targetfile[n-4], ".mnc") == 0)
    {
    ReadMINCImage(targetImage, targetMatrix, targetfile);
    }
  else
    {
    ReadDICOMImage(targetImage, targetMatrix, targetfile);
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
  sourceProperty->SetCheckerboardSpacing(40,40);

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

  renderWindow->SetSize(720,720);

  double bounds[6], center[4];
  targetImage->GetBounds(bounds);
  center[0] = 0.5*(bounds[0] + bounds[1]);
  center[1] = 0.5*(bounds[2] + bounds[3]);
  center[2] = 0.5*(bounds[4] + bounds[5]);
  center[3] = 1.0;
  targetMatrix->MultiplyPoint(center, center);

  vtkCamera *camera = renderer->GetActiveCamera();
  renderer->ResetCamera();
  camera->SetFocalPoint(center);
  camera->ParallelProjectionOn();
  camera->SetParallelScale(132);
  SetViewFromMatrix(renderer, istyle, targetMatrix);
  renderer->ResetCameraClippingRange();

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
  registration->SetInitializerTypeToCentered();
  registration->SetTransformTypeToRigid();
  //registration->SetTransformTypeToScaleTargetAxes();
  //registration->SetTransformTypeToAffine();
  //registration->SetMetricTypeToSquaredDifference();
  registration->SetMetricTypeToNormalizedMutualInformation();
  //registration->SetMetricTypeToNormalizedCrossCorrelation();
  //registration->SetMetricTypeToNeighborhoodCorrelation();
  registration->SetInterpolatorType(interpolatorType);
  registration->SetJointHistogramSize(numberOfBins,numberOfBins);
  registration->SetMetricTolerance(1e-4);
  registration->SetTransformTolerance(transformTolerance);
  registration->SetMaximumNumberOfIterations(500);

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

  for (;;)
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
    vtkSmartPointer<vtkMNITransformWriter> writer =
      vtkSmartPointer<vtkMNITransformWriter>::New();
    writer->SetFileName(xfmfile);
    writer->SetTransform(registration->GetTransform());
    writer->Update();
    }

  // -------------------------------------------------------
  // allow user to interact

  interactor->Start();

  return 1;
}
