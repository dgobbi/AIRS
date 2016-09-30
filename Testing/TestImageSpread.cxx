/*=========================================================================

  Module: ImageSpread.cxx

=========================================================================*/
// Test the ImageSpread class
//
// The command line arguments are:
// -I        => run in interactive mode


#include <vtkSmartPointer.h>
#include <vtkCamera.h>
#include <vtkImageData.h>
#include <vtkImageProperty.h>
#include <vtkImageReader2.h>
#include <vtkImageSlice.h>
#include <vtkImageSliceMapper.h>
#include <vtkROIStencilSource.h>
#include <vtkInteractorStyleImage.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkVersion.h>

#include <vtkTestUtilities.h>

#include "AIRSConfig.h"
#include "vtkImageSpread.h"

int main(int argc, char *argv[])
{
  vtkSmartPointer<vtkRenderWindowInteractor> iren =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  vtkSmartPointer<vtkInteractorStyleImage> style =
    vtkSmartPointer<vtkInteractorStyleImage>::New();
  style->SetInteractionModeToImageSlicing();
  vtkSmartPointer<vtkRenderWindow> renWin =
    vtkSmartPointer<vtkRenderWindow>::New();
  iren->SetRenderWindow(renWin);
  iren->SetInteractorStyle(style);

  char *fname =
    vtkTestUtilities::ExpandDataFileName(argc, argv, "Data/headsq/quarter");

  vtkSmartPointer<vtkImageReader2> reader =
    vtkSmartPointer<vtkImageReader2>::New();
  reader->SetDataByteOrderToLittleEndian();
  reader->SetDataExtent(0, 63, 0, 63, 2, 4);
  reader->SetDataSpacing(3.2, 3.2, 1.5);
  reader->SetFilePrefix(fname);

  delete [] fname;

  vtkSmartPointer<vtkROIStencilSource> stencil =
    vtkSmartPointer<vtkROIStencilSource>::New();
  stencil->SetOutputWholeExtent(0, 63, 0, 63, 2, 4);
  stencil->SetOutputOrigin(0.0, 0.0, 0.0);
  stencil->SetOutputSpacing(3.2, 3.2, 1.5);
  stencil->SetShapeToCylinderZ();
  stencil->SetBounds(0.0, 201.6, 50.0, 150.0, 0.0, 140.0);
  stencil->Update();

  for (int i = 0; i < 12; i++)
    {
    int j = i % 4;
    int k = i / 4;
    vtkSmartPointer<vtkRenderer> renderer =
      vtkSmartPointer<vtkRenderer>::New();
    vtkCamera *camera = renderer->GetActiveCamera();
    renderer->SetBackground(0.0, 0.0, 0.0);
    renderer->SetViewport(k/3.0, j/4.0, (k + 1)/3.0, (j + 1)/4.0);
    renWin->AddRenderer(renderer);

    vtkSmartPointer<vtkImageSpread> spread =
      vtkSmartPointer<vtkImageSpread>::New();
    spread->SetInputConnection(reader->GetOutputPort());
    spread->SetInputConnection(1, stencil->GetOutputPort());
    spread->SetNumberOfIterations(i);
    spread->Update();

    vtkSmartPointer<vtkImageSliceMapper> imageMapper =
      vtkSmartPointer<vtkImageSliceMapper>::New();
    imageMapper->SetInputConnection(spread->GetOutputPort());
    imageMapper->BorderOn();
    imageMapper->SliceFacesCameraOn();
    imageMapper->SliceAtFocalPointOn();

    double point[3] = { 100.8, 100.8, 5.25 };
    camera->SetFocalPoint(point);
    point[2] += 500.0;
    camera->SetPosition(point);
    camera->SetViewUp(0.0, 1.0, 0.0);
    camera->ParallelProjectionOn();
    camera->SetParallelScale(3.2*32);

    vtkSmartPointer<vtkImageSlice> image =
      vtkSmartPointer<vtkImageSlice>::New();
    image->SetMapper(imageMapper);
    image->GetProperty()->SetColorWindow(2047.0);
    image->GetProperty()->SetColorLevel(1023.5);
    renderer->AddViewProp(image);
    }

  renWin->SetSize(192, 256);

  iren->Initialize();
  renWin->Render();
  iren->Start();

  return EXIT_SUCCESS;
}
