/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageIslandRemoval.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageIslandRemoval - Select image regions by size
// .SECTION Description
// vtkImageIslandRemoval will identify connected regions in an image
// and will either keep or remove each region based on the size of that
// region.  A connected region is a group of voxels that are all within
// the supplied lower and upper intenstity thresholds and that furthermore
// are all connected to each other through their faces (or edges, in the
// case of 2D images).  The intensity thresholds are set similarly to
// vtkImageThreshold.  The scalar type of the output is the same as the
// input.
// .SECTION see also
// vtkImageThreshold vtkImageThresholdConnectivity
// .SECTION Thanks
// Thanks to David Gobbi at the Dept. of Radiology and and Dept. of Clinical
// Neurosciences, University of Calgary, for providing this class.

#ifndef vtkImageIslandRemoval_h
#define vtkImageIslandRemoval_h

#include "vtkImageSegmentationModule.h" // For export macro
#include "vtkImageAlgorithm.h"

class vtkPoints;
class vtkImageData;
class vtkImageStencilData;

class VTKIMAGESEGMENTATION_EXPORT vtkImageIslandRemoval :
  public vtkImageAlgorithm
{
public:
  static vtkImageIslandRemoval *New();
  vtkTypeMacro(vtkImageIslandRemoval, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Description:
  // Values greater than or equal to this threshold will be filled.
  void ThresholdByUpper(double thresh);

  // Description:
  // Values less than or equal to this threshold will be filled.
  void ThresholdByLower(double thresh);

  // Description:
  // Values within this range will be filled, where the range inludes
  // values that are exactly equal to the lower and upper thresholds.
  void ThresholdBetween(double lower, double upper);

  // Description:
  // Set the minimum island size to keep.  Size is in voxels.
  // Default value is zero.
  vtkSetMacro(SmallestIsland, vtkIdType);
  vtkGetMacro(SmallestIsland, vtkIdType);

  // Description:
  // Set the maximum island size to keep.  Size is in voxels.
  // Default value is VTK_ID_MAX.
  vtkSetMacro(LargestIsland, vtkIdType);
  vtkGetMacro(LargestIsland, vtkIdType);

  // Description:
  // Interpret LargestIsland and SmallestIsland as indices into a
  // list of all island sorted by size from largest to smallest.
  // So if LargestIsland is 0 and SmallestIsland is 1, then the
  // two largest islands will be kept.
  vtkSetMacro(IslandsSortedBySize, int);
  vtkBooleanMacro(IslandsSortedBySize, int);
  vtkGetMacro(IslandsSortedBySize, int);

  // Description:
  // Replace the filled region by the value set by SetInValue().
  vtkSetMacro(ReplaceIn, int);
  vtkGetMacro(ReplaceIn, int);
  vtkBooleanMacro(ReplaceIn, int);

  // Description:
  // If ReplaceIn is set, the filled region will be replaced by this value.
  void SetInValue(double val);
  vtkGetMacro(InValue, double);

  // Description:
  // Replace outside the filled region by the value set by SetOutValue().
  vtkSetMacro(ReplaceOut, int);
  vtkGetMacro(ReplaceOut, int);
  vtkBooleanMacro(ReplaceOut, int);

  // Description:
  // If ReplaceOut is set, outside the fill will be replaced by this value.
  void SetOutValue(double val);
  vtkGetMacro(OutValue, double);

  // Description:
  // Get the Upper and Lower thresholds.
  vtkGetMacro(UpperThreshold, double);
  vtkGetMacro(LowerThreshold, double);

  // Description:
  // Limit the flood to a range of slices in the specified direction.
  vtkSetVector2Macro(SliceRangeX, int);
  vtkGetVector2Macro(SliceRangeX, int);
  vtkSetVector2Macro(SliceRangeY, int);
  vtkGetVector2Macro(SliceRangeY, int);
  vtkSetVector2Macro(SliceRangeZ, int);
  vtkGetVector2Macro(SliceRangeZ, int);

  // Description:
  // Specify a stencil that will be used to limit the flood fill to
  // an arbitrarily-shaped region of the image.
  void SetStencilData(vtkImageStencilData *stencil);
  void SetStencil(vtkImageStencilData *stencil) {
    this->SetStencilData(stencil); }
  vtkImageStencilData *GetStencil();

  // Description:
  // For multi-component images, you can set which component will be
  // used for the threshold checks.
  vtkSetMacro(ActiveComponent,int);
  vtkGetMacro(ActiveComponent,int);

  // Description:
  // Override the MTime to account for the seed points.
  unsigned long GetMTime() override;

protected:
  vtkImageIslandRemoval();
  ~vtkImageIslandRemoval();

  double UpperThreshold;
  double LowerThreshold;
  double InValue;
  double OutValue;
  double IslandValue;
  int ReplaceIn;
  int ReplaceOut;
  int ReplaceIsland;

  vtkIdType LargestIsland;
  vtkIdType SmallestIsland;
  int IslandsSortedBySize;

  int SliceRangeX[2];
  int SliceRangeY[2];
  int SliceRangeZ[2];

  int ActiveComponent;

  vtkImageData *ImageMask;

  void ComputeInputUpdateExtent(int inExt[6], int outExt[6]);

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int RequestUpdateExtent(vtkInformation *, vtkInformationVector **,
                          vtkInformationVector *) override;
  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *) override;

private:
  vtkImageIslandRemoval(const vtkImageIslandRemoval&);  // Not implemented.
  void operator=(const vtkImageIslandRemoval&);  // Not implemented.
};

#endif
