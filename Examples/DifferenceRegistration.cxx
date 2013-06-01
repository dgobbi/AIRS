/*=========================================================================

Program:   Atamai Image Registration and Segmentation
Module:    DifferenceRegistration.cxx

   This software is distributed WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/

// This example registers two images (assumed to be of the same patient)
// and then subtracts first image from the second.

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

#include <vtkProperty.h>
#include <vtkDataSetMapper.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageSlice.h>
#include <vtkImageStack.h>
#include <vtkImageResliceMapper.h>
#include <vtkImageProperty.h>
#include <vtkImageHistogramStatistics.h>
#include <vtkImageMathematics.h>
#include <vtkImageShiftScale.h>
#include <vtkImageThreshold.h>
#include <vtkImageThresholdConnectivity.h>
#include <vtkImageContinuousDilate3D.h>
#include <vtkImageContinuousErode3D.h>

#include <vtkTimerLog.h>

#include <vtkImageRegistration.h>
#include <vtkImageMRIBrainExtractor.h>
#include <vtkImageIslandRemoval.h>
#include <vtkPolyDataToImageStencil.h>

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
#if 0
  const char *sourcefile2 = "";
  const char *targetfile2 = "";
#endif
  bool display = true;
  bool longitudinal = false;

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
  targetfile = argv[argi];
  sourcefile = argv[argi + 1];
#if 0
  if (argc > argi + 1)
    {
    targetfile2 = argv[argi + 2];
    sourcefile2 = argv[argi + 3];
    longitudinal = false;
    }
#endif

  // -------------------------------------------------------
  // parameters for registration

  int interpolatorType = vtkImageRegistration::Linear;
  double transformTolerance = 0.1; // tolerance on transformation result
  int numberOfBins = 64; // for Mattes' mutual information
  double initialBlurFactor = 4.0;

  // -------------------------------------------------------
  // load the images

  int n = 0;

  cerr << "source, target " << sourcefile << " " << targetfile << "\n";

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

#if 0
  // -------------------------------------------------------
  // Read the source image
  vtkSmartPointer<vtkImageData> sourceImage2 =
    vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkMatrix4x4> sourceMatrix2 =
    vtkSmartPointer<vtkMatrix4x4>::New();
  n = strlen(sourcefile2);
  if (n > 4 && strcmp(&sourcefile2[n-4], ".mnc") == 0)
    {
    ReadMINCImage(sourceImage2, sourceMatrix2, sourcefile2);
    }
  else if (n > 4)
    {
    ReadDICOMImage(sourceImage2, sourceMatrix2, sourcefile2);
    }

  // Read the target image
  vtkSmartPointer<vtkImageData> targetImage2 =
    vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkMatrix4x4> targetMatrix2 =
    vtkSmartPointer<vtkMatrix4x4>::New();
  n = strlen(targetfile2);
  if (n > 4 && strcmp(&targetfile[n-4], ".mnc") == 0)
    {
    ReadMINCImage(targetImage2, targetMatrix2, targetfile2);
    }
  else if (n > 4)
    {
    ReadDICOMImage(targetImage2, targetMatrix2, targetfile2);
    }
#endif

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
  //istyle->SetInteractionModeToImage3D();
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
  vtkSmartPointer<vtkImageHistogramStatistics> autoRange =
    vtkSmartPointer<vtkImageHistogramStatistics>::New();
  autoRange->SetInput(sourceImage);
  autoRange->Update();
  autoRange->GetAutoRange(sourceRange);

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
  autoRange->SetInput(targetImage);
  autoRange->Update();
  autoRange->GetAutoRange(targetRange);

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
  registration->SetMetricTypeToNormalizedMutualInformation();
  //registration->SetMetricTypeToNormalizedCrossCorrelation();
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
    //if (stage == 1)
      {
      blurFactor /= 2.0;
      if (blurFactor < 1.9)
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
  // do the subtraction
  double iScale = targetRange[1]/sourceRange[1]*1.2;
  //  (targetRange[1] - targetRange[0])/(sourceRange[1] - sourceRange[0]);
  double iShift = 0.0; //targetRange[0]/iScale - sourceRange[0];

  vtkSmartPointer<vtkImageSincInterpolator> sincInterpolator =
    vtkSmartPointer<vtkImageSincInterpolator>::New();
  sincInterpolator->SetWindowFunctionToBlackman();

  vtkSmartPointer<vtkImageReslice> resample =
    vtkSmartPointer<vtkImageReslice>::New();
  resample->SetInput(sourceImage);
  resample->SetInformationInput(targetImage);
  resample->SetInterpolator(sincInterpolator);
  resample->SetResliceTransform(registration->GetTransform()->GetInverse());
  resample->Update();

  vtkSmartPointer<vtkImageShiftScale> shiftScale =
    vtkSmartPointer<vtkImageShiftScale>::New();
  shiftScale->SetShift(iShift);
  shiftScale->SetScale(iScale);
  shiftScale->SetInput(resample->GetOutput());
  shiftScale->ClampOverflowOn();
  shiftScale->Update();

  vtkSmartPointer<vtkImageMathematics> difference =
    vtkSmartPointer<vtkImageMathematics>::New();
  difference->SetOperationToSubtract();
  difference->SetInput(0, shiftScale->GetOutput());
  difference->SetInput(1, targetImage);
  difference->Update();

  vtkSmartPointer<vtkImageContinuousDilate3D> dilate =
    vtkSmartPointer<vtkImageContinuousDilate3D>::New();
  dilate->SetKernelSize(3,3,3);
  dilate->SetInput(difference->GetOutput());
  dilate->Update();

/*
  vtkSmartPointer<vtkPoints> pts =
    vtkSmartPointer<vtkPoints>::New();
  vtkDataSet *pd = difference->GetOutput();
  vtkIdType np = pd->GetNumberOfPoints();
  pts->SetNumberOfPoints(np);
  for (vtkIdType ip = 0; ip < np; ip++)
    {
    double pt[3];
    pd->GetPoint(ip, pt);
    pts->SetPoint(ip, pt);
    }
*/

  vtkSmartPointer<vtkImageThreshold> thresh =
    vtkSmartPointer<vtkImageThreshold>::New();
  thresh->SetInput(dilate->GetOutput());
  thresh->ThresholdBetween(-5000, 2500.0);
  //thresh->SetSeedPoints(pts);
  thresh->ReplaceInOn();
  thresh->ReplaceOutOn();
  thresh->SetInValue(1);
  thresh->SetOutValue(0);
  thresh->Update();

  vtkSmartPointer<vtkImageMathematics> multiply =
    vtkSmartPointer<vtkImageMathematics>::New();
  multiply->SetOperationToMultiply();
  multiply->SetInput(0, difference->GetOutput());
  multiply->SetInput(1, thresh->GetOutput());
  multiply->Update();

  vtkSmartPointer<vtkImageSincInterpolator> imageBlurKernel =
    vtkSmartPointer<vtkImageSincInterpolator>::New();
  imageBlurKernel->SetWindowFunctionToBlackman();
  //imageBlurKernel->SetBlurFactors(5.0, 5.0, 5.0);

  vtkSmartPointer<vtkImageResize> imageBlur =
    vtkSmartPointer<vtkImageResize>::New();
  imageBlur->SetInput(difference->GetOutput());
  //imageBlur->SetInput(multiply->GetOutput());
  imageBlur->SetInterpolator(imageBlurKernel);
  imageBlur->InterpolateOn();
  imageBlur->Update();

  vtkSmartPointer<vtkImageContinuousDilate3D> dilate2 =
    vtkSmartPointer<vtkImageContinuousDilate3D>::New();
  dilate2->SetKernelSize(3,3,3);
  dilate2->SetInput(multiply->GetOutput());
  dilate2->Update();

  vtkSmartPointer<vtkImageContinuousErode3D> erode2 =
    vtkSmartPointer<vtkImageContinuousErode3D>::New();
  erode2->SetKernelSize(5,5,5);
  erode2->SetInput(dilate2->GetOutput());
  erode2->Update();

  vtkSmartPointer<vtkImageContinuousDilate3D> dilate3 =
    vtkSmartPointer<vtkImageContinuousDilate3D>::New();
  dilate3->SetKernelSize(3,3,3);
  dilate3->SetInput(erode2->GetOutput());
  dilate3->Update();

  double differenceRange[2];
  autoRange->SetInput(difference->GetOutput());
  autoRange->Update();
  autoRange->GetAutoRange(differenceRange);
  cout << differenceRange[0] << ", " << differenceRange[1] << "\n";
  differenceRange[0] = 0.0;

  vtkSmartPointer<vtkPoints> spoints =
    vtkSmartPointer<vtkPoints>::New();
  spoints->InsertNextPoint(88, 135, 80);

  vtkSmartPointer<vtkCellArray> verts =
    vtkSmartPointer<vtkCellArray>::New();
  verts->InsertNextCell(1);
  verts->InsertCellPoint(0);

  vtkSmartPointer<vtkPolyData> pdata =
    vtkSmartPointer<vtkPolyData>::New();
  pdata->SetPoints(spoints);
  pdata->SetVerts(verts);

  vtkSmartPointer<vtkDataSetMapper> pmapper =
    vtkSmartPointer<vtkDataSetMapper>::New();
  pmapper->SetInput(pdata);

  vtkSmartPointer<vtkActor> pactor =
    vtkSmartPointer<vtkActor>::New();
  pactor->SetMapper(pmapper);
  pactor->SetUserMatrix(targetMatrix);
  pactor->GetProperty()->SetRepresentationToPoints();

  renderer->AddViewProp(pactor);

  vtkSmartPointer<vtkImageIslandRemoval> seg =
    vtkSmartPointer<vtkImageIslandRemoval>::New();
  seg->SetInput(dilate3->GetOutput());
  //seg->SetSeedPoints(spoints);
  seg->ThresholdByUpper(0.10*differenceRange[1]);
  seg->ReplaceOutOn();
  seg->SetOutValue(0);
  //seg->IslandsSortedBySizeOn();
  //seg->SetSmallestIsland(1);
  //seg->SetLargestIsland(1);
  //seg->SetSmallestIsland(100000);
  //seg->SetLargestIsland(1000000);
  seg->Update();

  // -------------------------------------------------------
  // strip the skull
  vtkSmartPointer<vtkImageMRIBrainExtractor> bet =
    vtkSmartPointer<vtkImageMRIBrainExtractor>::New();
  bet->SetInput(targetImage);
  bet->Update();

  vtkSmartPointer<vtkPolyDataToImageStencil> maskMaker =
    vtkSmartPointer<vtkPolyDataToImageStencil>::New();
  maskMaker->SetInput(bet->GetBrainMesh());
  maskMaker->SetInformationInput(targetImage);
  maskMaker->Update();

  vtkSmartPointer<vtkImageReslice> imageMasker =
    vtkSmartPointer<vtkImageReslice>::New();
  //imageMasker->SetInput(difference->GetOutput());
  imageMasker->SetInput(imageBlur->GetOutput());
  //imageMasker->SetStencil(maskMaker->GetOutput());
  imageMasker->Update();

  // -------------------------------------------------------
  // display the result
  sourceMapper->SetInput(imageMasker->GetOutput());
  //sourceMapper->SetSlabThickness(100);
  //sourceMapper->SetSlabTypeToMax();

  sourceProperty->SetInterpolationTypeToLinear();
  sourceProperty->SetColorWindow((differenceRange[1]-differenceRange[0]));
  sourceProperty->SetColorLevel(0.5*(differenceRange[0]+differenceRange[1]));
  sourceProperty->CheckerboardOff();

  // -------------------------------------------------------
  // allow user to interact

  interactor->Start();

  return 1;
}
