/*=========================================================================

  Program:   Atamai Image Registration and Segmentation
  Module:    StrainGlyph.cxx

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

#include "vtkITKXFMReader.h"
#include "vtkMNITransformReader.h"
#include "vtkNIFTIReader.h"
#include "vtkNIFTIWriter.h"

#include <vtkSmartPointer.h>

#include <vtkMath.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageReslice.h>
#include <vtkImageSlice.h>
#include <vtkImageStack.h>
#include <vtkImageResliceMapper.h>
#include <vtkImageProperty.h>
#include <vtkDataSetMapper.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>
#include <vtkArrowSource.h>
#include <vtkLineSource.h>
#include <vtkTensorGlyph.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkGeneralTransform.h>
#include <vtkGridTransform.h>
#include <vtkTransform.h>
#include <vtkImageSincInterpolator.h>
#include <vtkErrorCode.h>
#include <vtkTimerLog.h>
#include <vtkStringArray.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkAppendPolyData.h>
#include <vtkLookupTable.h>
#include <vtkPNGWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkWindowToImageFilter.h>

void printUsage(const char *cmdname)
{
    cout << "Usage: " << cmdname << " image.nii tensors.nii"
         << endl;
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

    /*
    vtkSmartPointer<vtkImageGaussianInterpolator> interp =
      vtkSmartPointer<vtkImageGaussianInterpolator>::New();
    interp->SetBlurFactors(blurFactors);

    // smooth the grid with a Gaussian
    vtkSmartPointer<vtkImageResize> smooth =
      vtkSmartPointer<vtkImageResize>::New();
    smooth->SetInputConnection(reader->GetOutputPort());
    smooth->SetInterpolator(interp);
    smooth->Update();
    */
    vtkSmartPointer<vtkImageGaussianSmooth> smooth =
      vtkSmartPointer<vtkImageGaussianSmooth>::New();
    smooth->SetInputConnection(reader->GetOutputPort());
    smooth->SetRadiusFactors(4.5, 4.5, 4.5);
    smooth->SetStandardDeviation(
      0.399*blurFactors[0], 0.399*blurFactors[1], 0.399*blurFactors[2]);
    smooth->Update();

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

    vtkSmartPointer<vtkGridTransform> gt =
      vtkSmartPointer<vtkGridTransform>::New();
    // use linear to match ANTS?
    gt->SetInterpolationModeToCubic();
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

void ComputeGreensStrain(
  const double F[3][3], double G[3][3])
{
  vtkMath::Transpose3x3(F, G);
  vtkMath::Multiply3x3(G, F, G);

  G[0][0] -= 1.0;
  G[1][1] -= 1.0;
  G[2][2] -= 1.0;

  G[0][0] *= 0.5;
  G[0][1] *= 0.5;
  G[0][2] *= 0.5;
  G[1][0] *= 0.5;
  G[1][1] *= 0.5;
  G[1][2] *= 0.5;
  G[2][0] *= 0.5;
  G[2][1] *= 0.5;
  G[2][2] *= 0.5;
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

int main (int argc, char *argv[])
{
  if (argc < 3)
    {
    printUsage(argv[0]);
    return EXIT_FAILURE;
    }

  // -------------------------------------------------------
  // the files
  const char *imagefile;
  const char *tensorfile;

  int argi = 1;
  imagefile = argv[argi++];
  tensorfile = argv[argi++];

  vtkSmartPointer<vtkNIFTIReader> reader =
    vtkSmartPointer<vtkNIFTIReader>::New();
  reader->SetFileName(imagefile);
  reader->Update();
  vtkImageData *image = reader->GetOutput();

  vtkSmartPointer<vtkNIFTIReader> tensorReader =
    vtkSmartPointer<vtkNIFTIReader>::New();
  tensorReader->SetFileName(tensorfile);
  tensorReader->Update();
  vtkImageData *tensorImage = reader->GetOutput();

  // -------------------------------------------------------
  // read the affine transform to get the "average" tensor
  vtkSmartPointer<vtkGeneralTransform> xform =
    vtkSmartPointer<vtkGeneralTransform>::New();

  while (argi < argc)
    {
    static const double spacing[3] = { 1.0, 1.0, 1.0 }; //dummy
    const char *fname = argv[argi++];
    //size_t n = strlen(fname);
    //if (n > 4 && strcmp(fname[n-4], ".txt") == 0)
      {
      strain_read_transform(xform, fname, false, spacing);
      }
    }

  // -------------------------------------------------------
  // display the tensor glyphs

  double tensorSpacing[3], tensorOrigin[3];
  int tensorExtent[6];
  tensorImage->GetSpacing(tensorSpacing);
  tensorImage->GetOrigin(tensorOrigin);
  tensorImage->GetExtent(tensorExtent);
  int middleSlice = (tensorExtent[4] + tensorExtent[5])/2;
  int rfactor = 20;
  int trim = 1;
  tensorSpacing[0] = tensorSpacing[0]*rfactor;
  tensorSpacing[1] = tensorSpacing[1]*rfactor;
  tensorSpacing[2] = tensorSpacing[2]*1.0;
  tensorOrigin[0] += tensorSpacing[0]*trim;
  tensorOrigin[1] += tensorSpacing[1]*trim;
  tensorOrigin[2] += tensorSpacing[2]*middleSlice;
  tensorExtent[1] = tensorExtent[0] +
    (tensorExtent[1] - tensorExtent[0] + 1)/rfactor - 2*trim;
  tensorExtent[3] = tensorExtent[2] +
    (tensorExtent[3] - tensorExtent[2] + 1)/rfactor - 2*trim;
  tensorExtent[4] = 0;
  tensorExtent[5] = 0;

  vtkSmartPointer<vtkImageSincInterpolator> interp =
    vtkSmartPointer<vtkImageSincInterpolator>::New();
  interp->SetWindowFunctionToBlackman();
  interp->SetBlurFactors(rfactor, rfactor, 1.0);

  vtkSmartPointer<vtkImageReslice> tensorReslice =
    vtkSmartPointer<vtkImageReslice>::New();
  tensorReslice->SetInputConnection(tensorReader->GetOutputPort());
  tensorReslice->SetInterpolator(interp);
  tensorReslice->SetOutputOrigin(tensorOrigin);
  tensorReslice->SetOutputSpacing(tensorSpacing);
  tensorReslice->SetOutputExtent(tensorExtent);
  tensorReslice->Update();

  vtkSmartPointer<vtkImageData> tensors =
    vtkSmartPointer<vtkImageData>::New();
  tensors->CopyStructure(tensorReslice->GetOutput());
  tensors->GetPointData()->SetTensors(
    tensorReslice->GetOutput()->GetPointData()->GetScalars());

  //vtkSmartPointer<vtkSphereSource> glyphSource =
  //  vtkSmartPointer<vtkSphereSource>::New();
  //glyphSource->SetPhiResolution(19);
  //glyphSource->SetThetaResolution(9);

  vtkSmartPointer<vtkArrowSource> arrowSource =
    vtkSmartPointer<vtkArrowSource>::New();
  arrowSource->SetTipResolution(19);
  arrowSource->SetShaftResolution(19);
  arrowSource->SetTipRadius(0.03);
  arrowSource->SetShaftRadius(0.03);
  arrowSource->SetTipLength(0.90);

  vtkSmartPointer<vtkLineSource> lineSource =
    vtkSmartPointer<vtkLineSource>::New();
  lineSource->SetPoint1(0.0, 0.0, 0.0);
  lineSource->SetPoint2(1.0, 0.0, 0.0);

  vtkSmartPointer<vtkAppendPolyData> glyphSource =
    vtkSmartPointer<vtkAppendPolyData>::New();
  glyphSource->AddInputConnection(arrowSource->GetOutputPort());
  glyphSource->AddInputConnection(lineSource->GetOutputPort());

  double glyphScale = 2.0*rfactor;
  vtkSmartPointer<vtkTensorGlyph> glyph =
    vtkSmartPointer<vtkTensorGlyph>::New();
  glyph->SetInput(tensors);
  glyph->SetSource(glyphSource->GetOutput());
  glyph->SetScaleFactor(glyphScale);
  glyph->ThreeGlyphsOn();
  glyph->SymmetricOn();
  glyph->ColorGlyphsOn();
  glyph->SetColorModeToEigenvalues();

  vtkSmartPointer<vtkLookupTable> glyphTable =
    vtkSmartPointer<vtkLookupTable>::New();
  glyphTable->SetTableRange(-0.1, 0.1);
  double color[4] = { 0.0, 0.0, 0.0, 1.0 };
  for (int c = 0; c < 256; c++)
    {
    color[0] = (c <= 127 ? (127 - c)/127.0 : 0.0);
    color[1] = (c > 127 ? (c - 128)/127.0 : 0.0);
    glyphTable->SetTableValue(c, color);
    }

  vtkSmartPointer<vtkDataSetMapper> glyphMapper =
    vtkSmartPointer<vtkDataSetMapper>::New();
  glyphMapper->SetInputConnection(glyph->GetOutputPort());
  //glyphMapper->ScalarVisibilityOff();
  glyphMapper->SetLookupTable(glyphTable);
  glyphMapper->UseLookupTableScalarRangeOn();

  vtkSmartPointer<vtkActor> glyphActor =
    vtkSmartPointer<vtkActor>::New();
  glyphActor->SetMapper(glyphMapper);
  glyphActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
  glyphActor->GetProperty()->SetAmbient(0.5);
  glyphActor->GetProperty()->SetDiffuse(0.5);
  glyphActor->SetPosition(0.0, 0.0, -100.0);

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
  renderWindow->SetMultiSamples(16);

  vtkSmartPointer<vtkImageSlice> imageActor =
    vtkSmartPointer<vtkImageSlice>::New();
  vtkSmartPointer<vtkImageResliceMapper> imageMapper =
    vtkSmartPointer<vtkImageResliceMapper>::New();
  vtkSmartPointer<vtkImageProperty> imageProperty =
    vtkSmartPointer<vtkImageProperty>::New();

  imageMapper->SetInputConnection(reader->GetOutputPort());
  imageMapper->SliceAtFocalPointOn();
  imageMapper->SliceFacesCameraOn();
  imageMapper->ResampleToScreenPixelsOff();

  double sourceRange[2];
  image->GetScalarRange(sourceRange);
  imageProperty->SetInterpolationTypeToLinear();
  imageProperty->SetColorWindow((sourceRange[1]-sourceRange[0]));
  imageProperty->SetColorLevel(0.5*(sourceRange[0]+sourceRange[1]));

  imageActor->SetMapper(imageMapper);
  imageActor->SetProperty(imageProperty);

  vtkSmartPointer<vtkImageStack> imageStack =
    vtkSmartPointer<vtkImageStack>::New();
  imageStack->AddImage(imageActor);

  renderer->AddViewProp(imageStack);
  renderer->AddViewProp(glyphActor);
  renderer->SetBackground(0,0,0);

  // -------------------------------------------------------
  // make a tensor glyph from the transform
  if (argc > 3)
    {
    double point[3] = { 0.0, 0.0, 0.0 };
    double tensor[3][3];
    xform->InternalTransformDerivative(point, point, tensor);
    ComputeGreensStrain(tensor, tensor);

    vtkSmartPointer<vtkFloatArray> tensorArray =
      vtkSmartPointer<vtkFloatArray>::New();
    tensorArray->SetNumberOfComponents(9);
    tensorArray->InsertNextTuple(*tensor);

    cerr << tensor[0][0] << "\n";

    vtkSmartPointer<vtkPoints> tensorPoints =
      vtkSmartPointer<vtkPoints>::New();
    tensorPoints->InsertNextPoint(
      tensorExtent[1]*tensorSpacing[0] + tensorOrigin[0],
      tensorExtent[3]*tensorSpacing[1] + tensorOrigin[1],
      tensorExtent[5]*tensorSpacing[2] + tensorOrigin[2]);

    vtkIdType numPointIds = 1;
    static const vtkIdType pointIds[1] = { 0 };
    vtkSmartPointer<vtkCellArray> tensorVerts =
      vtkSmartPointer<vtkCellArray>::New();
    tensorVerts->InsertNextCell(numPointIds, pointIds);

    vtkSmartPointer<vtkPolyData> tensorData =
      vtkSmartPointer<vtkPolyData>::New();
    tensorData->SetPoints(tensorPoints);
    tensorData->SetVerts(tensorVerts);
    tensorData->GetPointData()->SetTensors(tensorArray);

    vtkSmartPointer<vtkTensorGlyph> glyph2 =
      vtkSmartPointer<vtkTensorGlyph>::New();
    glyph2->SetInput(tensorData);
    glyph2->SetSource(glyphSource->GetOutput());
    glyph2->SetScaleFactor(glyphScale);
    glyph2->ThreeGlyphsOn();
    glyph2->SymmetricOn();
    glyph2->ColorGlyphsOn();
    glyph2->SetColorModeToEigenvalues();

    vtkSmartPointer<vtkDataSetMapper> glyph2Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    glyph2Mapper->SetInputConnection(glyph2->GetOutputPort());
    //glyph2Mapper->ScalarVisibilityOff();
    glyphMapper->SetLookupTable(glyphTable);
    glyphMapper->UseLookupTableScalarRangeOn();

    vtkSmartPointer<vtkActor> glyph2Actor =
      vtkSmartPointer<vtkActor>::New();
    glyph2Actor->SetMapper(glyph2Mapper);
    glyph2Actor->GetProperty()->SetColor(0.0, 1.0, 0.0);
    glyph2Actor->GetProperty()->SetAmbient(0.5);
    glyph2Actor->GetProperty()->SetDiffuse(0.5);
    glyph2Actor->SetPosition(0.0, 0.0, -100.0);

    renderer->AddViewProp(glyph2Actor);
    }

  // -------------------------------------------------------
  renderWindow->SetSize(1024,1024);

  double bounds[6], center[3];
  image->GetBounds(bounds);
  center[0] = 0.5*(bounds[0] + bounds[1]);
  center[1] = 0.5*(bounds[2] + bounds[3]);
  center[2] = 0.5*(bounds[4] + bounds[5]);

  static const double viewRight[3] = { 1.0, 0.0, 0.0 };
  static const double viewUp[3] = { 0.0, -1.0, 0.0 };

  vtkCamera *camera = renderer->GetActiveCamera();
  renderer->ResetCamera();
  camera->SetFocalPoint(center);
  camera->ParallelProjectionOn();
  camera->SetParallelScale(0.5*(bounds[3] - bounds[2]));
  istyle->SetCurrentRenderer(renderer);
  istyle->SetImageOrientation(viewRight, viewUp);
  renderer->ResetCameraClippingRange();

  renderWindow->Render();

  // -------------------------------------------------------
  // allow user to interact

  WriteScreenshot(renderWindow, "tmp.png");

  interactor->Start();

  return 1;
}
