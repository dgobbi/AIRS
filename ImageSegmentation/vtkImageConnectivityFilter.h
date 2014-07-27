/*=========================================================================

  Program:   Atamai Image Registration and Segmentation
  Module:    vtkImageConnectivityFilter.h

  Copyright (c) 2014 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageConnectivityFilter - Label voxels by connectivity
// .SECTION Description
// vtkImageConnectivityFilter will identify and label connected regions
// in an image.  A connected region is a group of voxels in an image
// that have similar values and that are connected to each other through
// through their faces (or, for 2D images, pixels that are connected through
// their edges).  The labels for the regions can be applied based on the
// relative sizes of the regions, or can be explictly set based on the
// supplied seeds.
// .SECTION see also
// vtkImageThresholdConnectivity

#ifndef __vtkImageConnectivityFilter_h
#define __vtkImageConnectivityFilter_h

#include "vtkImageAlgorithm.h"

class vtkIdTypeArray;
class vtkIntArray;
class vtkPointSet;
class vtkImageData;
class vtkImageStencilData;

class VTK_EXPORT vtkImageConnectivityFilter :
  public vtkImageAlgorithm
{
public:
  static vtkImageConnectivityFilter *New();
  vtkTypeMacro(vtkImageConnectivityFilter, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Enum constants for SetLabelMode().
  enum LabelModeEnum {
    SeedScalar = 0,
    ConstantValue = 1,
    SizeRank = 2
  };

  // Description:
  // Enum constants for SetExtractionMode().
  enum ExtractionModeEnum {
    SeededRegions = 0,
    AllRegions = 1,
    LargestRegion = 2
  };

  // Description:
  // The input for seed locations (input port 1).
  // This must be a vtkPolyData or some other kind of vtkPointSet.
  // Each point in the supplied data set will be used as a seed.
  void SetSeedConnection(vtkAlgorithmOutput *port);
  vtkAlgorithmOutput *GetSeedConnection();
  void SetSeedData(vtkPointSet *data);

  // Description:
  // The input for a stencil (input port 2).
  // The labels will be restricted to the region inside the stencil,
  // as if no voxels existed outside the stencil.  This allows you to
  // apply this filter within an arbitrary region of interest.
  void SetStencilConnection(vtkAlgorithmOutput *port);
  vtkAlgorithmOutput *GetStencilConnection();
  void SetStencilData(vtkImageStencilData *data);

  // Description:
  // Set the scalar type for the output label image.
  // This should be one of UnsignedChar, Short, UnsignedShort, or Int
  // depending on how many labels are expected.  The default is UnsignedChar,
  // which allows for 255 label values.
  void SetLabelScalarTypeToUnsignedChar() {
    this->SetLabelScalarType(VTK_UNSIGNED_CHAR); }
  void SetLabelScalarTypeToShort() {
    this->SetLabelScalarType(VTK_SHORT); }
  void SetLabelScalarTypeToUnsignedShort() {
    this->SetLabelScalarType(VTK_UNSIGNED_SHORT); }
  void SetLabelScalarTypeToInt() {
    this->SetLabelScalarType(VTK_INT); }
  const char *GetLabelScalarTypeAsString();
  vtkSetMacro(LabelScalarType, int);
  vtkGetMacro(LabelScalarType, int);

  // Description:
  // Set the mode for applying labels to the output.
  // Labeling by SeedScalar uses the scalars of the seeds, or if
  // there are no scalars, then the regions will be labeled consecutively
  // starting at 1. Labeling by SizeRank means that the largest region is
  // labeled 1 and other regions are labeled consecutively in order of
  // decreasing size.  If there is a tie, then the seed point ID is used
  // as a tiebreaker.  Finally, Constant means that all regions will have
  // the value of SetLabelConstantValue().
  void SetLabelModeToSeedScalar() { this->SetLabelMode(SeedScalar); }
  void SetLabelModeToConstantValue() { this->SetLabelMode(ConstantValue); }
  void SetLabelModeToSizeRank() { this->SetLabelMode(SizeRank); }
  const char *GetLabelModeAsString();
  vtkSetMacro(LabelMode, int);
  vtkGetMacro(LabelMode, int);

  // Description:
  // Set which regions to output from the filter.
  void SetExtractionModeToSeededRegions(){
    this->SetExtractionMode(SeededRegions); }
  void SetExtractionModeToAllRegions() {
    this->SetExtractionMode(AllRegions); }
  void SetExtractionModeToLargestRegion() {
    this->SetExtractionMode(LargestRegion); }
  const char *GetExtractionModeAsString();
  vtkSetMacro(ExtractionMode, int);
  vtkGetMacro(ExtractionMode, int);

  // Description:
  // The label used when LabelMode is ConstantValue.
  vtkSetMacro(LabelConstantValue, int);
  vtkGetMacro(LabelConstantValue, int);

  // Description:
  // Get the number of extracted regions.
  vtkIdType GetNumberOfExtractedRegions();

  // Description:
  // Get the label used for each extracted region.
  vtkIntArray *GetExtractedRegionLabels() {
    return this->ExtractedRegionLabels; }

  // Desciption:
  // Get the size of each extracted region.
  vtkIdTypeArray *GetExtractedRegionSizes() {
    return this->ExtractedRegionSizes; }

  // Description:
  // Get the PointId of the seed for each region.
  // If no seed was used, the PointId will be -1.
  vtkIdTypeArray *GetExtractedRegionSeedIds() {
    return this->ExtractedRegionSeedIds; }

  // Description:
  // Set the size range for the extracted regions.
  // Only regions that have sizes within the specified range will be present
  // in the output.  The default range is (1, VTK_ID_MAX).
  vtkSetVector2Macro(SizeRange, vtkIdType);
  vtkGetVector2Macro(SizeRange, vtkIdType);

  // Description:
  // Set the scalar range used to define potential regions.
  // Only voxels with values that are within this range will be considered
  // for region membership.  This is an inclusive range, meaning that the
  // upper and lower limits are considered to be within the range.  The
  // default range goes from 0.5 to VTK_DOUBLE_MAX.
  vtkSetVector2Macro(ScalarRange, double);
  vtkGetVector2Macro(ScalarRange, double);

  // Description:
  // For multi-component input images, select which component to use.
  vtkSetMacro(ActiveComponent, int);
  vtkGetMacro(ActiveComponent, int);

  // Description:
  // Limit the operation to a slice range in each direction.
  // This is most often used to limit the connectivity operation to
  // just a single slice of a 3D image.
  vtkSetVector2Macro(SliceRangeX, int);
  vtkGetVector2Macro(SliceRangeX, int);
  vtkSetVector2Macro(SliceRangeY, int);
  vtkGetVector2Macro(SliceRangeY, int);
  vtkSetVector2Macro(SliceRangeZ, int);
  vtkGetVector2Macro(SliceRangeZ, int);

protected:
  vtkImageConnectivityFilter();
  ~vtkImageConnectivityFilter();

  int LabelMode;
  int ExtractionMode;

  double ScalarRange[2];
  vtkIdType SizeRange[2];
  int LabelConstantValue;
  int ActiveComponent;
  int LabelScalarType;

  int SliceRangeX[2];
  int SliceRangeY[2];
  int SliceRangeZ[2];

  vtkIntArray *ExtractedRegionLabels;
  vtkIdTypeArray *ExtractedRegionSizes;
  vtkIdTypeArray *ExtractedRegionSeedIds;

  void ComputeInputUpdateExtent(int inExt[6], int outExt[6]);

  virtual int FillInputPortInformation(int port, vtkInformation *info);
  virtual int RequestInformation(
    vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int RequestUpdateExtent(
    vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int RequestData(
    vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  vtkImageConnectivityFilter(const vtkImageConnectivityFilter&);  // Not implemented.
  void operator=(const vtkImageConnectivityFilter&);  // Not implemented.
};

#endif
