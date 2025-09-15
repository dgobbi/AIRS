/*=========================================================================

  Program:   Visualization Toolkit
  Module:    ImageConnectivityFilter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// Test the ImageConnectivityFilter class
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
#include "vtkImageConnectivityFilter.h"

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

  vtkSmartPointer<vtkPoints> seedPoints =
    vtkSmartPointer<vtkPoints>::New();
  seedPoints->InsertNextPoint(1, 1, 5.25);
  seedPoints->InsertNextPoint(100.8, 100.8, 5.25);
  vtkSmartPointer<vtkUnsignedCharArray> seedScalars =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  seedScalars->InsertNextValue(2);
  seedScalars->InsertNextValue(5);
  vtkSmartPointer<vtkPolyData> seedData =
    vtkSmartPointer<vtkPolyData>::New();
  seedData->SetPoints(seedPoints);
  seedData->GetPointData()->SetScalars(seedScalars);

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

    vtkSmartPointer<vtkImageConnectivityFilter> connectivity =
      vtkSmartPointer<vtkImageConnectivityFilter>::New();
    connectivity->SetInputConnection(reader->GetOutputPort());
    //connectivity->SetSeedData(seedData);
    //connectivity->SetExtractionModeToSeededRegions();
    connectivity->SetLabelModeToSizeRank();
    if (k == 0)
    {
      connectivity->SetScalarRange(0, 800);
    }
    else if (k == 1)
    {
      connectivity->SetScalarRange(1200, 4095);
    }
    else
    {
      connectivity->SetScalarRange(800, 1200);
    }
    connectivity->GenerateRegionExtentsOn();

    // test a previous bug where OutputExtent != InputExtent cause a crash.
    int extent[6] = { 0, 63, 0, 63, 3, 3 };
    connectivity->UpdateExtent(extent);

    vtkIdTypeArray *sizeArray = connectivity->GetExtractedRegionSizes();
    vtkIdTypeArray *idArray = connectivity->GetExtractedRegionSeedIds();
    vtkIdTypeArray *labelArray = connectivity->GetExtractedRegionLabels();
    vtkIntArray *extentArray = connectivity->GetExtractedRegionExtents();
    vtkIdType rn = connectivity->GetNumberOfExtractedRegions();
    cout << "info";
    for (vtkIdType r = 0; r < rn; r++)
    {
      cout << " (" << idArray->GetValue(r) << ","
           << labelArray->GetValue(r) << ","
           << sizeArray->GetValue(r) << ",["
           << extentArray->GetValue(6*r) << ","
           << extentArray->GetValue(6*r+1) << ","
           << extentArray->GetValue(6*r+2) << ","
           << extentArray->GetValue(6*r+3) << ","
           << extentArray->GetValue(6*r+4) << ","
           << extentArray->GetValue(6*r+5) << "])";
    }
    cout << "\n";

    vtkSmartPointer<vtkImageSliceMapper> imageMapper =
      vtkSmartPointer<vtkImageSliceMapper>::New();
    imageMapper->SetInputConnection(connectivity->GetOutputPort());
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
    image->GetProperty()->SetColorWindow(6);
    image->GetProperty()->SetColorLevel(3);
    renderer->AddViewProp(image);
  }

  renWin->SetSize(192, 256);

  iren->Initialize();
  renWin->Render();
  iren->Start();

  return EXIT_SUCCESS;
}
