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
#include <vtkImageSincInterpolator.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkMath.h>

#include <vtkMINCImageReader.h>
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

  vtkImageData *image = reader->GetOutput();

  // get the data
  data->CopyStructure(image);
  data->GetPointData()->PassData(image->GetPointData());

  // generate the matrix
  matrix->DeepCopy(reader->GetDirectionCosines());
}

};

int main (int argc, char *argv[])
{
  if (argc < 3)
    {
    cout << "Usage 1: " << argv[0] << " source.mnc target.mnc" << endl;
    cout << "Usage 2: " << argv[0] << " dicomdir1/ dicomdir2/" << endl;
    return EXIT_FAILURE;
    }

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
  n = strlen(argv[1]);
  if (n > 4 && strcmp(&argv[1][n-4], ".mnc") == 0)
    {
    ReadMINCImage(sourceImage, sourceMatrix, argv[1]);
    }
  else
    {
    ReadDICOMImage(sourceImage, sourceMatrix, argv[1]);
    }

  // Read the target image
  vtkSmartPointer<vtkImageData> targetImage =
    vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkMatrix4x4> targetMatrix =
    vtkSmartPointer<vtkMatrix4x4>::New();
  n = strlen(argv[2]);
  if (n > 4 && strcmp(&argv[2][n-4], ".mnc") == 0)
    {
    ReadMINCImage(targetImage, targetMatrix, argv[2]);
    }
  else
    {
    ReadDICOMImage(targetImage, targetMatrix, argv[2]);
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

  istyle->SetInteractionModeToImage3D();
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
  camera->SetPosition(center[0], center[1], center[2] - 100);
  camera->SetViewUp(0.0, 1.0, 0.0);
  camera->ParallelProjectionOn();
  camera->SetParallelScale(132);
  renderer->ResetCameraClippingRange();

  renderWindow->Render();

  // -------------------------------------------------------
  // prepare for registration

  // display the images
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

  // blur source image with Hamming-windowed sinc
  vtkSmartPointer<vtkImageSincInterpolator> sourceBlurKernel =
    vtkSmartPointer<vtkImageSincInterpolator>::New();
  sourceBlurKernel->SetWindowFunctionToHamming();

  // reduce the source resolution
  vtkSmartPointer<vtkImageReslice> sourceBlur =
    vtkSmartPointer<vtkImageReslice>::New();
  sourceBlur->SetInput(sourceImage);
  sourceBlur->SetInterpolator(sourceBlurKernel);

  // blur target with Hamming-windowed sinc
  vtkSmartPointer<vtkImageSincInterpolator> targetBlurKernel =
    vtkSmartPointer<vtkImageSincInterpolator>::New();
  targetBlurKernel->SetWindowFunctionToHamming();

  // keep target at full resolution
  vtkSmartPointer<vtkImageReslice> targetBlur =
    vtkSmartPointer<vtkImageReslice>::New();
  targetBlur->SetInput(targetImage);
  targetBlur->SetInterpolator(targetBlurKernel);

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
  //registration->SetTransformTypeToSimilarity();
  //registration->SetMetricTypeToNormalizedCrossCorrelation();
  registration->SetMetricTypeToNormalizedMutualInformation();
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

  // do a multi-resolution registration
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
      registration->SetTransformTolerance(minSpacing);
      }
    else
      {
      registration->SetInterpolatorType(interpolatorType);
      registration->SetTransformTolerance(transformTolerance);
      }
    if (blurFactor < 1.1)
      {
      // full resolution: no blurring or resampling
      sourceBlur->SetInterpolator(0);
      sourceBlur->SetInterpolationModeToNearestNeighbor();
      sourceBlur->SetOutputSpacing(sourceSpacing);
      sourceBlur->Update();

      targetBlur->SetInterpolator(0);
      targetBlur->SetInterpolationModeToNearestNeighbor();
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
      interactor->Render();
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
  // allow user to interact

  interactor->Start();

  return 1;
}
