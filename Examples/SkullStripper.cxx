/** -*- c++ -*- 
 *
 * \file   skullStripper.cxx
 * \date   Mon Apr 20 09:30:45 2015
 *
 * \copyright 
 * Copyright (c) 2015 Liangfu Chen <liangfu.chen@nlpr.ia.ac.cn>.
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms are permitted
 * provided that the above copyright notice and this paragraph are
 * duplicated in all such forms and that any documentation,
 * advertising materials, and other materials related to such
 * distribution and use acknowledge that the software was developed
 * by the Brainnetome Center & NLPR at Institute of Automation, CAS. The 
 * name of the Brainnetome Center & NLPR at Institute of Automation, CAS 
 * may not be used to endorse or promote products derived
 * from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * \brief  skull stripping utility
 */
 
#include "vtkSmartPointer.h"
#include "vtkNIFTIReader.h"
#include "vtkNIFTIWriter.h"
#include "vtkImageMRIBrainExtractor.h"
#include "vtkImageData.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyData.h"

int main(int argc, char * argv[])
{
  char filename[1024];
  if (argc<2){fprintf(stderr,"require args.\n");return -1;}else {strcpy(filename,argv[1]);}
  
  vtkSmartPointer<vtkNIFTIReader> reader = vtkSmartPointer<vtkNIFTIReader>::New();
  reader->SetFileName(filename);
  reader->Update();

  vtkImageData * imageData = reader->GetOutput();

  vtkSmartPointer<vtkImageMRIBrainExtractor> stripper=vtkSmartPointer<vtkImageMRIBrainExtractor>::New();
  stripper->SetInputData( imageData );

  // A kludge to handle low axial resolution images
  double spacing[3];
  imageData->GetSpacing(spacing);
  float zres = spacing[2];

  float bt;
  if (zres > 1.5){
	bt = 0.50;
  }else{
	bt = 0.70;
  }

  // Optimized parameters
  stripper->SetRMin(8.0);
  stripper->SetRMax(10.0);
  stripper->SetD1(7.0);
  stripper->SetD2(3.0);
  stripper->SetBT(bt);
  stripper->Update();

  // Get the basename
  // char basename[1024];
  // sprintf(basename,"output000");
  std::string basename = "output000";

  // Write the mesh
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetFileName((basename + "_mesh.vtk").c_str());
  writer->SetInputData( stripper->GetBrainMesh() );
  writer->Write();

  // Write the segemented image
  vtkSmartPointer<vtkNIFTIWriter> writer2 = vtkSmartPointer<vtkNIFTIWriter>::New();
  writer2->SetInputData( stripper->GetOutput() );
  writer2->SetFileName((basename + "_mask.nii").c_str());
  writer2->Write();

  return 0;
}
