/*=========================================================================

  Program:   Atamai Classes for VTK
  Module:    vtkFrameFinder.cxx

Copyright (c) 2005 Atamai, Inc.
All rights reserved.

Use, modification and redistribution of the software, in source or
binary forms, are permitted provided that the following terms and
conditions are met:

1) Redistribution of the source code, in verbatim or modified
   form, must retain the above copyright notice, this license,
   the following disclaimer, and any notices that refer to this
   license and/or the following disclaimer.

2) Redistribution in binary form must include the above copyright
   notice, a copy of this license and the following disclaimer
   in the documentation or with other materials provided with the
   distribution.

3) Modified copies of the source code must be clearly marked as such,
   and must not be misrepresented as verbatim copies of the source code.

THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE SOFTWARE "AS IS"
WITHOUT EXPRESSED OR IMPLIED WARRANTY INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE.  IN NO EVENT SHALL ANY COPYRIGHT HOLDER OR OTHER PARTY WHO MAY
MODIFY AND/OR REDISTRIBUTE THE SOFTWARE UNDER THE TERMS OF THIS LICENSE
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, LOSS OF DATA OR DATA BECOMING INACCURATE
OR LOSS OF PROFIT OR BUSINESS INTERRUPTION) ARISING IN ANY WAY OUT OF
THE USE OR INABILITY TO USE THE SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES.

=========================================================================*/

#include "vtkFrameFinder.h"

// VTK header files
#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkImageData.h"
#include "vtkPolyData.h"
#include "vtkMatrix4x4.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkVersion.h"
#include "vtkMath.h"
#include "vtkSmartPointer.h"

#include "vtkTemplateAliasMacro.h"
// turn off 64-bit ints when templating over all types
# undef VTK_USE_INT64
# define VTK_USE_INT64 0
# undef VTK_USE_UINT64
# define VTK_USE_UINT64 0

#include <math.h>
#include <string.h>

#include <stack>
#include <map>
#include <vector>
#include <utility>
#include <limits>

vtkStandardNewMacro(vtkFrameFinder);
vtkCxxSetObjectMacro(vtkFrameFinder,DICOMPatientMatrix,vtkMatrix4x4);

//----------------------------------------------------------------------------
vtkFrameFinder::vtkFrameFinder()
{
  this->ImageToFrameMatrix = vtkMatrix4x4::New();
  this->DICOMPatientMatrix = 0;

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkFrameFinder::~vtkFrameFinder()
{
  if (this->DICOMPatientMatrix)
    {
    this->DICOMPatientMatrix->Delete();
    }
  if (this->ImageToFrameMatrix)
    {
    this->ImageToFrameMatrix->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkFrameFinder::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "ImageToFrameMatrix: "
     << this->ImageToFrameMatrix << "\n";

  os << indent << "DICOMPatientMatrix: "
     << this->DICOMPatientMatrix << "\n";
}

//----------------------------------------------------------------------------
int vtkFrameFinder::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
  return 1;
}

//----------------------------------------------------------------------------
int vtkFrameFinder::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

//----------------------------------------------------------------------------
void vtkFrameFinder::SetInput(vtkDataObject* input)
{
#if VTK_MAJOR_VERSION <= 5
  if (input)
    {
    this->SetInputConnection(0, input->GetProducerPort());
    }
  else
    {
    // Setting a NULL input removes the connection.
    this->SetInputConnection(0, 0);
    }
#else
  this->SetInputDataInternal(0, input);
#endif
}

//----------------------------------------------------------------------------
vtkDataObject* vtkFrameFinder::GetInput()
{
  return this->GetExecutive()->GetInputData(0, 0);
}

//----------------------------------------------------------------------------
vtkPolyData *vtkFrameFinder::GetOutput()
{
  return vtkPolyData::SafeDownCast(this->GetOutputDataObject(0));
}

//----------------------------------------------------------------------------
int vtkFrameFinder::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *vtkNotUsed(outputVector))
{
  return 1;
}

//----------------------------------------------------------------------------
int vtkFrameFinder::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector))
{
  int inExt[6];

  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), inExt);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inExt, 6);

  return 1;
}

//----------------------------------------------------------------------------
int vtkFrameFinder::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkImageData *input = vtkImageData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (!this->DICOMPatientMatrix)
    {
    vtkErrorMacro("A DICOMPatientMatrix must be given");
    return 0;
    }

  double xdir[4] = { 1.0, 0.0, 0.0, 0.0 };
  double ydir[4] = { 0.0, 1.0, 0.0, 0.0 };
  double zdir[4] = { 0.0, 0.0, 1.0, 0.0 };
  this->DICOMPatientMatrix->MultiplyPoint(xdir, xdir);
  this->DICOMPatientMatrix->MultiplyPoint(ydir, ydir);
  this->DICOMPatientMatrix->MultiplyPoint(zdir, zdir);

  if (fabs(xdir[0]) < 0.9 || fabs(ydir[1]) < 0.9)
    {
    vtkErrorMacro("Frame images must be axial");
    return 0;
    }

  /*
  // generate the frame data
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

  output->SetPoints(points);
  output->SetLines(cells);
  */

  // indicate if x or y is flipped
  double direction[3];
  direction[0] = xdir[0]; // DICOM +x and Leksell +x go to the right
  direction[1] = -ydir[1]; // DICOM +y is posterior, Leksell +y is anterior
  direction[2] = -zdir[2]; // DICOM +z is superior, Leksell +z is inferior

  // find the frame in the input
  this->FindFrame(input, output, direction,
    this->ImageToFrameMatrix);

  return 1;
}

//----------------------------------------------------------------------------
int vtkFrameFinder::ProcessRequest(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // generate the data oject
  if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_OBJECT()))
    {
    return 1;
    }
  // generate the data
  if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
    {
    return this->RequestData(request, inputVector, outputVector);
    }

  // execute information
  if (request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
    {
    return this->RequestInformation(request, inputVector, outputVector);
    }

  // propagate update extent
  if (request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
    {
    return this->RequestUpdateExtent(request, inputVector, outputVector);
    }

  return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

namespace {
//----------------------------------------------------------------------------
// Frame-finding code
//----------------------------------------------------------------------------

// Store a single pixel and its value
struct Pixel
{
  int x;
  int y;
  int z;
  int val;
};

// Store a centroid computed from one or more pixels
struct Blob
{
  int count;
  int slice;
  double x;
  double y;
  double val;
};

// Store a point
struct Point
{
  double x;
  double y;
  double z;
};

//----------------------------------------------------------------------------
// Functions for finding blobs of neighboring bright pixels

// Get a pixel, then zero it to remove it from future consideration
template <class T>
void GetAndZeroPixel(
  T* array, int x, int y, int z,
  int pixelsPerRow, int pixelsPerSlice, T *value)
{
  *value = array[z*pixelsPerSlice + y*pixelsPerRow + x];
  array[z*pixelsPerSlice + y*pixelsPerRow + x] = 0;
}

// Add pixel to basket if in bounds and above threshold
template <class T>
void AddPixelIfAboveThreshold(
  T* array, int x, int y, int z,
  int dataSizeX, int dataSizeY, int dataSizeZ,
  T lowerThresh, T upperThresh,
  struct Blob *theBlob, std::stack<Pixel> *basket)
{
  T value;

  if (x >= 0 && y >= 0 && z >= 0 &&
      x < dataSizeX && y < dataSizeY && z < dataSizeZ)
    {
    GetAndZeroPixel(array, x, y, z, dataSizeX, dataSizeY*dataSizeX, &value);
    if (value > lowerThresh)
      {
      if (value > upperThresh) { value = upperThresh; }
      Pixel pD = { x, y, z, static_cast<int>(value) };
      (*basket).push(pD);
      theBlob->count += 1;
      theBlob->x += value*x;
      theBlob->y += value*y;
      theBlob->val += value;
      }
    }
}

// Perform a flood-fill to search for a blob of bright pixels
template <class T>
Blob MakeBlob(
  T* array, int x, int y, int z,
  int dataSizeX, int dataSizeY, int dataSizeZ,
  int blobMaxSize, T lowerThresh, T upperThresh)
{
  struct Blob theBlob;

  theBlob.count = 0;
  theBlob.slice = z;
  theBlob.x = 0.0;
  theBlob.y = 0.0;
  theBlob.val = 0.0;
  Pixel curPt;
  std::stack<Pixel> basket;

  AddPixelIfAboveThreshold(array, x, y, z, dataSizeX, dataSizeY, dataSizeZ,
                           lowerThresh, upperThresh, &theBlob, &basket);

  while (basket.size() > 0)
    {
    curPt = basket.top();
    basket.pop();
    AddPixelIfAboveThreshold(array, curPt.x + 1, curPt.y, z,
                             dataSizeX, dataSizeY, dataSizeZ,
                             lowerThresh, upperThresh, &theBlob, &basket);
    AddPixelIfAboveThreshold(array, curPt.x + 1, curPt.y + 1, z,
                             dataSizeX, dataSizeY, dataSizeZ,
                             lowerThresh, upperThresh, &theBlob, &basket);
    AddPixelIfAboveThreshold(array, curPt.x, curPt.y + 1, z,
                             dataSizeX, dataSizeY, dataSizeZ,
                             lowerThresh, upperThresh, &theBlob, &basket);
    AddPixelIfAboveThreshold(array, curPt.x - 1, curPt.y + 1, z,
                             dataSizeX, dataSizeY, dataSizeZ,
                             lowerThresh, upperThresh, &theBlob, &basket);
    }

  if (theBlob.count > 1 && theBlob.count <= blobMaxSize)
    {
    theBlob.x /= theBlob.val;
    theBlob.y /= theBlob.val;
    theBlob.val /= theBlob.count;
    }
  else
    {
    theBlob.count = 0;
    }

  return theBlob;
}

// Scan an image volume for blobs
template <class T>
void ScanVolumeForBlobs(
  T* inPtr, const int dataSize[3],
  double lower, double upper, int blobMaxSize,
  std::vector<Blob> *returnList)
{
  size_t rowSize = dataSize[0];
  size_t sliceSize = rowSize*dataSize[1];
  size_t volumeSize = sliceSize*dataSize[2];

  // compute a threshold and find the blobs
  double minVal = static_cast<double>(std::numeric_limits<T>::min());
  double maxVal = static_cast<double>(std::numeric_limits<T>::max());
  if (lower < minVal) { lower = minVal; }
  if (lower > maxVal) { lower = maxVal; }
  T lowerThresh = static_cast<T>(lower);
  if (upper < minVal) { upper = minVal; }
  if (upper > maxVal) { upper = maxVal; }
  T upperThresh = static_cast<T>(upper);

  T *array = new T [volumeSize];
  memcpy(array, inPtr, volumeSize*sizeof(T));

  T *pixelPtr = array;
  for (int z = 0; z < dataSize[2]; z++)
    {
    for (int y = 0; y < dataSize[1]; y++)
      {
      for (int x = 0; x < dataSize[0]; x++)
        {
        if (*pixelPtr++ > lowerThresh)
          {
          Blob blob = MakeBlob(array, x, y, z, dataSize[0], dataSize[1],
                               dataSize[2], blobMaxSize,
                               lowerThresh, upperThresh);
          if (blob.count != 0)
            {
            returnList->push_back(blob);
            }
          }
        }
      }
    }

  delete [] array;
}

// Find the threshold that includes the given fraction of pixel values
template <class T>
void ComputePercentile(
  T* inPtr, const int dataSize[3], const int boxSize[3],
  double fraction, double *threshold)
{
  size_t rowSize = dataSize[0];
  size_t sliceSize = rowSize*dataSize[1];

  // compute the min/max values in the box
  size_t offset = ((dataSize[0] - boxSize[0])/2 +
                   rowSize*(dataSize[1] - boxSize[1])/2 +
                   sliceSize*(dataSize[2] - boxSize[2])/2);

  T *boxPixelPtr = inPtr + offset;
  T minVal = std::numeric_limits<T>::max();
  T maxVal = std::numeric_limits<T>::min();

  for (int z = 0; z < boxSize[2]; z++)
    {
    for (int y = 0; y < boxSize[1]; y++)
      {
      for (int x = 0; x < boxSize[0]; x++)
        {
        T v = *boxPixelPtr;
        minVal = (minVal < v ? minVal : v);
        maxVal = (maxVal > v ? maxVal : v);
        boxPixelPtr++;
        }
      boxPixelPtr += dataSize[0] - boxSize[0];
      }
    boxPixelPtr += rowSize*(dataSize[1] - boxSize[1]);
    }

  if (minVal == maxVal)
    {
    *threshold = static_cast<int>(maxVal);
    return;
    }

  // compute the histogram in the box
  const int histSize = 1024;
  int *histo = new int [histSize];
  double shift = -minVal;
  double scale = static_cast<double>(histSize - 1)/(maxVal - minVal);
  if (scale > 1.0) { scale = 1.0; }

  for (int i = 0; i < histSize; i++)
    {
    histo[i] = 0;
    }

  size_t total = 0;
  boxPixelPtr = inPtr + offset;
  for (int z = 0; z < boxSize[2]; z++)
    {
    for (int y = 0; y < boxSize[1]; y++)
      {
      for (int x = 0; x < boxSize[0]; x++)
        {
        T v = *boxPixelPtr;
        int bin = static_cast<int>((v + shift)*scale);
        histo[bin]++;
        total++;
        boxPixelPtr++;
        }
      boxPixelPtr += dataSize[0] - boxSize[0];
      }
    boxPixelPtr += rowSize*(dataSize[1] - boxSize[1]);
    }

  // exclude pixels above the fraction
  size_t fractionTotal = static_cast<size_t>(fraction*total);
  size_t sum = 0;
  T thresh = 0;
  for (int i = 0; i < histSize; i++)
    {
    sum += histo[i];
    thresh = (sum > fractionTotal ? thresh : i);
    histo[i] = 0;
    }

  *threshold = thresh/scale - shift;
  delete [] histo;
}

// An adapter to call ScanVolumeForBlobs on vtkImageData
void UpdateBlobs(
  vtkImageData *input, std::vector<Blob> *blobs)
{
  // get basic image information
  void *address = input->GetScalarPointer();
  int dataType = input->GetScalarType();
  double spacing[3];
  input->GetSpacing(spacing);
  int extent[6];
  input->GetExtent(extent);

  // maximum blob size to find is 24 square millimetres
  int maxBlobSize = static_cast<int>(24.0/(spacing[0]*spacing[1]));

  // size of image volume
  int dataSize[3];
  dataSize[0] = extent[1] - extent[0] + 1;
  dataSize[1] = extent[3] - extent[2] + 1;
  dataSize[2] = extent[5] - extent[4] + 1;

  // use a centred 20x20x10cm box to compute image thresholds
  int boxSize[3];
  boxSize[0] = static_cast<int>(200.0/spacing[0]);
  boxSize[1] = static_cast<int>(200.0/spacing[1]);
  boxSize[2] = static_cast<int>(100.0/spacing[2]);

  // limit box size
  if (boxSize[0] > dataSize[0]) { boxSize[0] = dataSize[0]; }
  if (boxSize[1] > dataSize[1]) { boxSize[1] = dataSize[1]; }
  if (boxSize[2] > dataSize[2]) { boxSize[2] = dataSize[2]; }

  // threshold at half of the 98th percentile within the box
  const double fthresh = 0.5;
  const double fraction = 0.98;
  double thresh98; // to be computed by ComputePercentile

  // find the 98th percentile threshold
  switch (dataType)
    {
    vtkTemplateAliasMacro(
      ComputePercentile(static_cast<VTK_TT *>(address),
                        dataSize, boxSize, fraction, &thresh98));
    }

  // call the blob-finding function
  switch (dataType)
    {
    vtkTemplateAliasMacro(
      ScanVolumeForBlobs(static_cast<VTK_TT *>(address), dataSize,
                         fthresh*thresh98, thresh98,
                         maxBlobSize, blobs));
    }
}

//----------------------------------------------------------------------------

//
class Histogram
{
public:
  Histogram(std::vector<Blob> *blobs, const double direction[3]);
  double GetMinimum() { return this->Minimum; }
  double GetMaximum() { return this->Maximum; }
  double GetAverage() { return this->Average; }
  bool CollapseTwo(double distance, double threshold, int maxWidth);
  bool CollapseOne(double position, double threshold, int maxWidth);
  bool CollapseOne(double threshold, int maxWidth);

  struct Bin
  {
    int minSlice;
    int maxSlice;
    double sum;
  };

  struct Cluster
  {
    int lowest;
    int highest;
    int minSlice;
    int maxSlice;
    double sum;
  };

  std::vector<Cluster> *GetClusters() { return &this->Clusters; }

private:
  void ReturnCluster(double keepThreshold, Cluster *cluster);
  void BuildClusters(double threshold, int maxWidth);

  std::map<int, Bin> Bins;
  std::vector<Cluster> Clusters;
  double Minimum;
  double Maximum;
  double Average;
};

// collapse the blob points into bins by x or y value
Histogram::Histogram(std::vector<Blob> *blobs, const double vec[3])
{
  std::map<int, Bin> *bins = &this->Bins;

  double sum = 0.0;
  double minval = 1e30;
  double maxval = -1e30;

  for (std::vector<Blob>::iterator it = blobs->begin();
       it != blobs->end();
       ++it)
    {
    // dot product with the direction vector, then round to int
    double x = vec[0]*it->x + vec[1]*it->y + vec[2]*it->slice;
    int j = static_cast<int>(x > 0 ? x + 0.5 : x - 0.5);
    std::map<int, Bin>::iterator jt = bins->find(j);
    if (jt == bins->end())
      {
      Bin bin = { it->slice, it->slice, it->val };
      bins->insert(std::make_pair(j, bin));
      sum += it->val;
      minval = (minval < it->val ? minval : it->val);
      maxval = (maxval > it->val ? maxval : it->val);
      }
    else
      {
      jt->second.minSlice = std::min(jt->second.minSlice, it->slice);
      jt->second.maxSlice = std::max(jt->second.minSlice, it->slice);
      jt->second.sum += it->val;
      sum += it->val;
      minval = (minval < it->val ? minval : it->val);
      maxval = (maxval > jt->second.sum ? maxval : jt->second.sum);
      }
    }

  this->Minimum = minval;
  this->Maximum = maxval;
  this->Average = sum/bins->size();
}

void Histogram::ReturnCluster(double threshold, Cluster *cluster)
{
  // Picks largest peak in histogram, returns tuple consisting
  // of range of peak and area of peak.
  // Removes bins from Bins
  // Multiple calls return the multiple largest peaks in the
  // histogram

  // first find the peak bin
  std::map<int, Bin> *bins = &this->Bins;
  int minSlice = 0;
  int maxSlice = 0;
  int maxCluster = 0;
  double maxClusterValue = 0;
  for (std::map<int, Bin>::iterator it = bins->begin();
       it != bins->end();
       ++it)
    {
    int i = it->first;
    double v = it->second.sum;
    if (v > maxClusterValue)
      {
      minSlice = it->second.minSlice;
      maxSlice = it->second.maxSlice;
      maxCluster = i;
      maxClusterValue= v;
      }
    }
  double sum = maxClusterValue;
  bins->erase(maxCluster);

  // then count down from the middle
  int lowest = maxCluster;
  for (;;)
    {
    std::map<int, Bin>::iterator it = bins->find(lowest - 1);
    if (it == bins->end()) { break; }
    double v = it->second.sum;
    if (v <= threshold) { break; }
    sum += v;
    minSlice = std::min(minSlice, it->second.minSlice);
    maxSlice = std::max(maxSlice, it->second.maxSlice);
    bins->erase(it);
    lowest--;
    }

  // then count up from the middle
  int highest = maxCluster;
  for (;;)
    {
    std::map<int, Bin>::iterator it = bins->find(highest + 1);
    if (it == bins->end()) { break; }
    double v = it->second.sum;
    if (v <= threshold) { break; }
    sum += v;
    minSlice = std::min(minSlice, it->second.minSlice);
    maxSlice = std::max(maxSlice, it->second.maxSlice);
    bins->erase(it);
    highest++;
    }

  cluster->lowest = lowest;
  cluster->highest = highest;
  cluster->minSlice = minSlice;
  cluster->maxSlice = maxSlice;
  cluster->sum = sum;
}

void Histogram::BuildClusters(double threshold, int maxWidth)
{
  std::vector<Cluster> *clusters = &this->Clusters;
  std::map<int, Bin> *bins = &this->Bins;

  clusters->clear();
  while (!bins->empty())
    {
    Cluster cluster;
    this->ReturnCluster(threshold, &cluster);
    if (cluster.highest - cluster.lowest <= maxWidth &&
        //cluster.maxSlice - cluster.minSlice >= 50 &&
        cluster.sum > threshold)
      {
      cerr << cluster.highest << " - " << cluster.lowest << " <= " << maxWidth << " : " << cluster.sum << " > " << 0.5*this->Average << " : " << cluster.maxSlice << " - " << cluster.minSlice << "\n";
      clusters->push_back(cluster);
      }
    }
}

bool Histogram::CollapseTwo(
  double distance, double threshold, int maxWidth)
{
  // distance = ideal distance between histogram peaks in
  // data coordinates

  this->BuildClusters(threshold, maxWidth);
  std::vector<Cluster> *clusters = &this->Clusters;

  double smallestDiff = 1e30;

  typedef std::vector<Cluster>::iterator ClusterIterator;
  std::vector<std::pair<ClusterIterator, ClusterIterator> > candidates;
  size_t bestCandidate = 0;

  for (ClusterIterator it = clusters->begin(); it != clusters->end(); ++it)
    {
    double binPos = 0.5*(it->lowest + it->highest);
    for (ClusterIterator jt = clusters->begin(); jt != it; ++jt)
      {
      double binPos2 = 0.5*(jt->lowest + jt->highest);
      double curDist = fabs(binPos - binPos2);
      double curDiff = fabs(curDist - distance);
      if (curDiff < static_cast<double>(maxWidth))
        {
        if (curDiff < smallestDiff)
          {
          smallestDiff = curDiff;
          bestCandidate = candidates.size();
          }
        if (it->lowest < jt->lowest)
          {
          candidates.push_back(std::make_pair(it, jt));
          }
        else
          {
          candidates.push_back(std::make_pair(jt, it));
          }
        }
      }
    }

  if (candidates.size() == 0)
    {
    return false;
    }

  // merge clusters that are within maxwidth of best clusters
  for (size_t i = 0; i < candidates.size(); i++)
    {
    if (i != bestCandidate)
      {
      ClusterIterator a = candidates[bestCandidate].first;
      ClusterIterator b = candidates[i].first;
      for (int j = 0; j < 2; j++)
        {
        int dist = a->highest - b->lowest;
        if (dist < 0) { dist = b->highest - a->lowest; }
        if (a != b && dist <= maxWidth)
          {
          cerr << "joining " << b->lowest << " " << b->highest << "\n";
          a->lowest = std::min(a->lowest, b->lowest);
          a->highest = std::max(a->highest, b->highest);
          a->minSlice = std::min(a->minSlice, b->minSlice);
          a->maxSlice = std::max(a->maxSlice, b->maxSlice);
          a->sum += b->sum;
          }
        a = candidates[bestCandidate].second;
        b = candidates[i].second;
        }
      }
    }

  Cluster cluster1 = *(candidates[bestCandidate].first);
  Cluster cluster2 = *(candidates[bestCandidate].second);

  cerr << "smallestDiff " << smallestDiff << "\n";
  cerr << "pos " << (cluster1.lowest + cluster1.highest)/2
       << " " << (cluster2.lowest + cluster2.highest)/2
       << " sum " << cluster1.sum << " " << cluster2.sum
       << " range1 " << cluster1.minSlice << " " << cluster1.maxSlice
       << " range2 " << cluster2.minSlice << " " << cluster2.maxSlice
       << "\n";

  clusters->clear();
  clusters->push_back(cluster1);
  clusters->push_back(cluster2);

  return true;
}

bool Histogram::CollapseOne(
  double position, double threshold, int maxWidth)
{
  this->BuildClusters(threshold, maxWidth);
  std::vector<Cluster> *clusters = &this->Clusters;
  cerr << "CollapseOne clusters " << clusters->size() << "\n";

  double smallestDiff = 1e30;

  typedef std::vector<Cluster>::iterator ClusterIterator;
  std::vector<ClusterIterator> candidates;
  size_t bestCandidate = 0;

  for (ClusterIterator it = clusters->begin(); it != clusters->end(); ++it)
    {
    double binPos = 0.5*(it->lowest + it->highest);
    double curDiff = fabs(binPos - position);
    cerr << "CollapseOne curDiff " << curDiff << " " << it->sum << "\n";
    if (curDiff < static_cast<double>(maxWidth))
      {
      if (curDiff < smallestDiff)
        {
        smallestDiff = curDiff;
        bestCandidate = candidates.size();
        }
      candidates.push_back(it);
      }
    }

  if (candidates.size() == 0)
    {
    return false;
    }

  // merge clusters that are within maxwidth of best cluster
  for (size_t i = 0; i < candidates.size(); i++)
    {
    if (i != bestCandidate)
      {
      ClusterIterator a = candidates[bestCandidate];
      ClusterIterator b = candidates[i];
      int dist = a->highest - b->lowest;
      if (dist < 0) { dist = b->highest - a->lowest; }
      if (a != b && dist <= maxWidth)
        {
        a->lowest = std::min(a->lowest, b->lowest);
        a->highest = std::max(a->highest, b->highest);
        a->minSlice = std::min(a->minSlice, b->minSlice);
        a->maxSlice = std::max(a->maxSlice, b->maxSlice);
        a->sum += b->sum;
        }
      }
    }

  Cluster cluster = *(candidates[bestCandidate]);

  cerr << "smallestDiff " << smallestDiff << "\n";
  cerr << "pos " << (cluster.lowest + cluster.highest)/2
       << " sum " << cluster.sum
       << " range " << cluster.minSlice << " " << cluster.maxSlice
       << "\n";

  clusters->clear();
  clusters->push_back(cluster);

  return true;
}

bool Histogram::CollapseOne(double threshold, int maxWidth)
{
  this->BuildClusters(threshold, maxWidth);
  std::vector<Cluster> *clusters = &this->Clusters;

  typedef std::vector<Cluster>::iterator ClusterIterator;
  double maxSum = -1e30;
  ClusterIterator bestCluster = clusters->end();

  for (ClusterIterator it = clusters->begin(); it != clusters->end(); ++it)
    {
    if (it->sum > maxSum)
      {
      maxSum = it->sum;
      bestCluster = it;
      }
    }

  if (bestCluster == clusters->end())
    {
    return false;
    }

  Cluster cluster = *bestCluster;

  cerr << "pos " << (cluster.lowest + cluster.highest)/2
       << " sum " << cluster.sum
       << " range " << cluster.minSlice << " " << cluster.maxSlice
       << "\n";

  clusters->clear();
  clusters->push_back(cluster);

  return true;
}


class FiducialBar
{
public:
  FiducialBar() {
    this->Points = 0;
    this->Eigenvector[0] = this->Eigenvector[1] = this->Eigenvector[2] = 0;
    this->CentreOfMass[0] = this->CentreOfMass[1] = this->CentreOfMass[2] = 0;
    this->Eigenvalue = 0; }

  void SetPoints(std::vector<Point> *p) {
    this->Points = p;
    this->CentreOfMass[0] = 0;
    this->CentreOfMass[1] = 0;
    this->CentreOfMass[2] = 0; }

  std::vector<Point> *GetPoints() {
    return this->Points; }

  void ClipFiducialEnds(double zMin, double zMax);

  bool ExtractLinesFromPoints();

  void CullOutliers(double maxDist);

  double GetLine(double p[3], double v[3]);

private:
  void ComputeCentreOfMass(double com[3]);

  std::vector<Point> *Points;
  double Eigenvalue;
  double Eigenvector[3];
  double CentreOfMass[3];
};

double FiducialBar::GetLine(double p[3], double v[3])
{
  p[0] = this->CentreOfMass[0];
  p[1] = this->CentreOfMass[1];
  p[2] = this->CentreOfMass[2];

  v[0] = this->Eigenvector[0];
  v[1] = this->Eigenvector[1];
  v[2] = this->Eigenvector[2];

  if (v[2] < 0)
    {
    v[0] = -v[0];
    v[1] = -v[1];
    v[2] = -v[2];
    }

  return fabs(this->Eigenvalue);
}

void FiducialBar::ClipFiducialEnds(double minZ, double maxZ)
{
  size_t j = 0;
  for (size_t i = 0; i < this->Points->size(); i++)
    {
    Point &pi = (*this->Points)[i];
    if (pi.z > minZ && pi.z < maxZ)
      {
      (*this->Points)[j++] = pi;
      }
    }
  this->Points->resize(j);
}

void FiducialBar::ComputeCentreOfMass(double com[3])
{
  double count = 0;
  com[0] = 0.0;
  com[1] = 0.0;
  com[2] = 0.0;

  for (std::vector<Point>::iterator it = this->Points->begin();
       it < this->Points->end();
       ++it)
    {
    com[0] += it->x;
    com[1] += it->y;
    com[2] += it->z;
    count++;
    }

  if (count != 0)
    {
    com[0] /= count;
    com[1] /= count;
    com[2] /= count;
    }
}

bool FiducialBar::ExtractLinesFromPoints()
{
  if (this->Points->size() == 0)
    {
    this->CentreOfMass[0] = 0.0;
    this->CentreOfMass[1] = 0.0;
    this->CentreOfMass[2] = 0.0;
    this->Eigenvalue = 0.0;
    this->Eigenvector[0] = 0.0;
    this->Eigenvector[1] = 0.0;
    this->Eigenvector[2] = 0.0;
    return false;
    }

  double C[3][3] = { {0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0} };
  this->ComputeCentreOfMass(this->CentreOfMass);
  double *Q = this->CentreOfMass;
  for (std::vector<Point>::iterator it = this->Points->begin();
       it < this->Points->end();
       ++it)
    {
    C[0][0] += (it->x - Q[0])*(it->x - Q[0]);
    C[0][1] += (it->x - Q[0])*(it->y - Q[1]);
    C[0][2] += (it->x - Q[0])*(it->z - Q[2]);
    C[1][0] += (it->y - Q[1])*(it->x - Q[0]);
    C[1][1] += (it->y - Q[1])*(it->y - Q[1]);
    C[1][2] += (it->y - Q[1])*(it->z - Q[2]);
    C[2][0] += (it->z - Q[2])*(it->x - Q[0]);
    C[2][1] += (it->z - Q[2])*(it->y - Q[1]);
    C[2][2] += (it->z - Q[2])*(it->z - Q[2]);
    }

  double eigenvalues[3];
  double eigenvectors[3][3];

  // use matrix format that can be used with Jacobi routine
  double *CA[3], *EA[3];
  CA[0] = C[0];
  CA[1] = C[1];
  CA[2] = C[2];
  EA[0] = eigenvectors[0];
  EA[1] = eigenvectors[1];
  EA[2] = eigenvectors[2];
  vtkMath::Jacobi(CA, eigenvalues, EA);

  this->Eigenvalue = eigenvalues[0];
  this->Eigenvector[0] = eigenvectors[0][0];
  this->Eigenvector[1] = eigenvectors[1][0];
  this->Eigenvector[2] = eigenvectors[2][0];

  double *v = this->Eigenvector;
  cerr << "(" << Q[0] << ", " << Q[1] << ", " << Q[2] << ") ";
  cerr << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")\n";

  return true;
}

void FiducialBar::CullOutliers(double maxDist)
{
  // Check each point in the bar for distance from the extracted line
  // If the distance exceeds maxDist, remove the point from the list.

  double maxDistSquared = maxDist*maxDist;
  double sum = 0.0;
  double *com = this->CentreOfMass;

  size_t j = 0;
  for (size_t i = 0; i < this->Points->size(); i++)
    {
    Point &pi = (*this->Points)[i];

    // compute distance from point to line with pythagoras
    double p[3];
    double *v = this->Eigenvector;
    p[0] = pi.x - com[0];
    p[1] = pi.y - com[1];
    p[2] = pi.z - com[2];
    double dotpp = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
    double dotpv = v[0]*p[0] + v[1]*p[1] + v[2]*p[2];
    double distSquared = dotpp - dotpv*dotpv;

    sum += distSquared;

    if (distSquared < maxDistSquared)
      {
      (*this->Points)[j++] = pi;
      }
    }

  // compute what the RMS was originally
  if (this->Points->size())
    {
    sum = sqrt(sum/static_cast<double>(this->Points->size()));
    }

  cerr << "RMS " << sum << " maxDist " << maxDist << "\n";

  if (sum > maxDist)
    {
    // RMS was too high, throw away all the points
    this->Points->clear();
    }
  else
    {
    // points too far from line have been discarded
    this->Points->resize(j);
    }
}

void LineIntersection(
  const double p1[3], const double v1[3],
  const double p2[3], const double v2[3],
  double p[3])
{
  double v3[3];
  v3[0] = p1[0] - p2[0];
  v3[1] = p1[1] - p2[1];
  v3[2] = p1[2] - p2[2];

  double a = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
  double b = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
  double c = v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2];
  double d = v1[0]*v3[0] + v1[1]*v3[1] + v1[2]*v3[2];
  double e = v2[0]*v3[0] + v2[1]*v3[1] + v2[2]*v3[2];
  double D = a*c - b*b;
  double t1 = 0.0;
  double t2 = 0.0;
  if (D > 1e-6)
    {
    t1 = (b*e - c*d)/D;
    t2 = (a*e - b*d)/D;
    }

  p[0] = 0.5*((p1[0] + v1[0]*t1) + (p2[0] + v2[0]*t2));
  p[1] = 0.5*((p1[1] + v1[1]*t1) + (p2[1] + v2[1]*t2));
  p[2] = 0.5*((p1[2] + v1[2]*t1) + (p2[2] + v2[2]*t2));
}

class FiducialPlate
{
public:
  FiducialPlate() {};

  bool LocateBars(
    std::vector<Blob> *blobs, const double hvec[3], const double dvec[3],
    const double origin[3], const double spacing[3],
    double barSeparation, double clusterThreshold, double clusterWidth,
    double zMin, double zMax);

  FiducialBar *GetBars() { return this->Bars; }

private:
  // diagonal bar first, then vertical bars
  FiducialBar Bars[3];
  std::vector<Point> BarPoints[3];
};

bool FiducialPlate::LocateBars(
  std::vector<Blob> *blobs, const double hvec[3], const double dvec[3],
  const double origin[3], const double spacing[3],
  double barSeparation, double clusterThreshold, double clusterWidth,
  double zLow, double zHigh)
{
  // which direction are we going in?
  double hspacing = spacing[1];
  if (fabs(hvec[0]) > fabs(hvec[1]))
    {
    hspacing = spacing[0];
    }

  // diagonal bar is longer, so use larger cluster width
  int clusterWidthD = static_cast<int>(clusterWidth/hspacing*1.4 + 0.5);
  int clusterWidthH = static_cast<int>(clusterWidth/hspacing + 0.5);

  // generate histograms along the horiz and diagonal directions
  Histogram hist1(blobs, dvec);
  Histogram hist2(blobs, hvec);

  // the cluster identification threshold
  clusterThreshold *= hist2.GetMaximum();

  // locate the diagonal bar
  if (!hist1.CollapseOne(clusterThreshold, clusterWidthD))
    {
    cerr << "diagonal failed\n";
    return false;
    }
  cerr << "diagonal success\n";

  // locate the two vertical bars
  if (!hist2.CollapseTwo(barSeparation, clusterThreshold, clusterWidthH))
    {
    cerr << "vertical failed\n";
    return false;
    }
  cerr << "vertical success\n";

  // get the clusters for the three bars
  std::vector<Histogram::Cluster>::iterator clusters[3];
  clusters[0] = hist1.GetClusters()->begin();
  clusters[1] = hist2.GetClusters()->begin();
  clusters[2] = hist2.GetClusters()->begin() + 1;

  // the histogram directions for the bars
  const double *vecs[3];
  vecs[0] = dvec;
  vecs[1] = hvec;
  vecs[2] = hvec;

  // find out which bar each blob belongs to
  for (std::vector<Blob>::iterator it = blobs->begin();
       it != blobs->end();
       ++it)
    {
    // convert blob position to 3D data coordinates
    Point p;
    p.x = origin[0] + spacing[0]*it->x;
    p.y = origin[1] + spacing[1]*it->y;
    p.z = origin[2] + spacing[2]*it->slice;

    // loop through the three bars (diagonal first!)
    for (int i = 0; i < 3; i++)
      {
      // compute the histogram index for the point
      const double *vec = vecs[i];
      double x = it->x*vec[0] + it->y*vec[1] + it->slice*vec[2];
      int idx = static_cast<int>(x > 0 ? x + 0.5 : x - 0.5);

      if (idx >= clusters[i]->lowest &&
          idx <= clusters[i]->highest)
        {
        this->BarPoints[i].push_back(p);
        break;
        }
      }
    }

  // do the computations for each bar
  bool result = true;
  for (int i = 0; i < 3; i++)
    {
    FiducialBar *bar = &this->Bars[i];
    bar->SetPoints(&this->BarPoints[i]);
    bar->ClipFiducialEnds(zLow, zHigh);
    bar->ExtractLinesFromPoints();
    bar->CullOutliers(2.0);
    if (!bar->ExtractLinesFromPoints())
      {
      cerr << "failed on bar " << i << "\n";
      result = false;
      }
    }

  return result;
}

void BuildMatrix(
  const double xvec[3], const double yvec[3], const double zvec[3],
  const double centre[3], const double frameCentre[3],
  const double direction[3], double matrix[16])
{
  double M[3][3];
  for (int j = 0; j < 3; j++)
    {
    M[0][j] = xvec[j];
    M[1][j] = yvec[j];
    M[2][j] = zvec[j];
    }

  vtkMath::Normalize(M[0]);
  vtkMath::Normalize(M[1]);
  vtkMath::Normalize(M[2]);

  vtkMath::Orthogonalize3x3(M, M);

  double dsign[3];
  dsign[0] = (direction[0] < 0 ? -1.0 : +1.0);
  dsign[1] = (direction[1] < 0 ? -1.0 : +1.0);
  dsign[2] = (direction[2] < 0 ? -1.0 : +1.0);

  for (int k = 0; k < 3; k++)
    {
    matrix[0 + k] = M[0][k]*dsign[0];
    matrix[4 + k] = M[1][k]*dsign[1];
    matrix[8 + k] = M[2][k]*dsign[2];
    matrix[12 + k] = 0.0;
    }

  double tc[3];
  tc[0] = matrix[0]*centre[0] + matrix[1]*centre[1] + matrix[2]*centre[2];
  tc[1] = matrix[4]*centre[0] + matrix[5]*centre[1] + matrix[6]*centre[2];
  tc[2] = matrix[8]*centre[0] + matrix[9]*centre[1] + matrix[10]*centre[2];

  matrix[3] = frameCentre[0] - tc[0];
  matrix[7] = frameCentre[1] - tc[1];
  matrix[11] = frameCentre[2] - tc[2];
  matrix[15] = 1.0;
}

bool PositionFrame(
  std::vector<Blob> *blobs, std::vector<Point> *points,
  const int extent[6], const double origin[3], const double spacing[3],
  const double direction[3], double matrix[16])
{
  cerr << "spacing " << spacing[0] << " " << spacing[1] << " " << spacing[2] << "\n";

  double plateSeparationX = 196.0;
  double plateSeparationY = 235.0;
  double barSeparation = 120.0;
  double plateClusterThreshold = 0.1;
  double barClusterThreshold = 0.05;
  double clusterWidth = 10.0;

  int clusterWidthX = static_cast<int>(clusterWidth/spacing[0] + 0.5);
  int clusterWidthY = static_cast<int>(clusterWidth/spacing[1] + 0.5);

  double xvec[3] = { 1.0, 0.0, 0.0 };
  double yvec[3] = { 0.0, 1.0, 0.0 };
  double zvec[3] = { 0.0, 0.0, 1.0 };
  double frameCentre[3] = { 100.0, 100.0, 100.0 };

  double centre[3];
  centre[0] = 0.5*(extent[0] + extent[1])*spacing[0] + origin[0];
  centre[1] = 0.5*(extent[2] + extent[3])*spacing[1] + origin[1];
  centre[2] = 0.5*(extent[4] + extent[5])*spacing[2] + origin[2];

  // create an initial matrix
  BuildMatrix(xvec, yvec, zvec, centre, frameCentre, direction, matrix);

  // find the side plates of the leksell frame
  Histogram xHistogram(blobs, xvec);
  if (!xHistogram.CollapseTwo(plateSeparationX/spacing[0],
                              plateClusterThreshold*xHistogram.GetAverage(),
                              clusterWidthX))
    {
    cerr << "could not find sides\n";
    return false;
    }
  std::vector<Histogram::Cluster> *xClusters = xHistogram.GetClusters();

  // for computing min and max slice for side plates
  int zMin = 10000;
  int zMax = -1;

  // for storing the blobs that belong to each of the side plates
  std::vector<Blob> plateBlobs[4];

  // find the blobs that are within the side-plate clusters,
  // or within plus or minus one pixel of these clusters
  for (std::vector<Blob>::iterator it = blobs->begin();
       it != blobs->end();
       ++it)
    {
    for (int j = 0; j < 2; j++)
      {
      double x = it->x;
      int xIdx = static_cast<int>(x > 0 ? x + 0.5 : x - 0.5);
      if (xIdx >= (*xClusters)[j].lowest &&
          xIdx <= (*xClusters)[j].highest)
        {
        zMin = (zMin < it->slice ? zMin : it->slice);
        zMax = (zMax > it->slice ? zMax : it->slice);
        plateBlobs[j].push_back(*it);
        break;
        }
      }
    }

  // we will be chopping 10% from the top and bottom of the plates
  // in order to better capture each fiducial bar in isolation from
  // other bars
  double zLow = (zMin + (zMax - zMin)*0.1)*spacing[2] + origin[2];
  double zHigh = (zMax - (zMax - zMin)*0.1)*spacing[2] + origin[2];

  // direction perpendicular to the diagonal bars
  double dvec[3];
  dvec[0] = 0.0;
  dvec[1] = sqrt(0.5);
  dvec[2] = -sqrt(0.5)*spacing[2]/spacing[1]*direction[1]*direction[2];

  // find the bars for each of the side plates
  FiducialPlate plates[4];
  bool foundPlate[4];
  for (int j = 0; j < 2; j++)
    {
    foundPlate[j] = plates[j].LocateBars(
      &plateBlobs[j], yvec, dvec, origin, spacing, barSeparation/spacing[1],
      barClusterThreshold, clusterWidth, zLow, zHigh);
    }

  if (!foundPlate[0] || !foundPlate[1])
    {
    cerr << "could not find the side plates\n";
    return false;
    }

  // these will be used to find the front/back plates
  double ycentre = 0.0;

  // will be computing the intersection points at the corners of the
  // plates and the average direction of the vertical bars
  double corners[8][3];
  zvec[0] = 0.0;
  zvec[1] = 0.0;
  zvec[2] = 0.0;

  for (int j = 0; j < 2; j++)
    {
    double barVectors[3][3];
    double barCentroids[3][3];
    double barWeights[3];

    for (int k = 0; k < 3; k++)
      {
      FiducialBar *bar = &plates[j].GetBars()[k];
      std::vector<Point> *barPoints = bar->GetPoints();
      points->insert(points->end(), barPoints->begin(), barPoints->end());
      barWeights[k] = bar->GetLine(barCentroids[k], barVectors[k]);

      // for the vertical bars only
      if (k > 0)
        {
        zvec[0] += barVectors[k][0]*barWeights[k];
        zvec[1] += barVectors[k][1]*barWeights[k];
        zvec[2] += barVectors[k][2]*barWeights[k];
        }
      }

    // compute points where diagonal bar intersects vertical bars
    LineIntersection(
      barCentroids[0], barVectors[0], barCentroids[1], barVectors[1],
      corners[2*j + 0]);
    LineIntersection(
      barCentroids[0], barVectors[0], barCentroids[2], barVectors[2],
      corners[2*j + 1]);

    // compute approximate centre in the anterior-posterior direction
    ycentre += 0.25*(barCentroids[1][1] + 0.5*barSeparation);
    ycentre += 0.25*(barCentroids[2][1] - 0.5*barSeparation);
    }

  // find the blobs within which to search for front/back plates
  double yPlatePos[2];
  yPlatePos[0] = ycentre - 0.5*plateSeparationY;
  yPlatePos[1] = ycentre + 0.5*plateSeparationY;
  double yPlateRange[2][2];
  yPlateRange[0][0] = (yPlatePos[0]-0.5*clusterWidth - origin[1])/spacing[1];
  yPlateRange[0][1] = (yPlatePos[0]+0.5*clusterWidth - origin[1])/spacing[1];
  yPlateRange[1][0] = (yPlatePos[1]-0.5*clusterWidth - origin[1])/spacing[1];
  yPlateRange[1][1] = (yPlatePos[1]+0.5*clusterWidth - origin[1])/spacing[1];

  cerr << "y plate positions " << yPlatePos[0] << ", " << yPlatePos[1] << "\n";

  std::vector<Blob> yBlobs[2];
  for (std::vector<Blob>::iterator it = blobs->begin();
       it != blobs->end();
       ++it)
    {
    for (int j = 0; j < 2; j++)
      {
      if (it->y > yPlateRange[j][0] && it->y < yPlateRange[j][1])
        {
        yBlobs[j].push_back(*it);
        break;
        }
      }
    }

  dvec[0] = sqrt(0.5);
  dvec[1] = 0.0;
  dvec[2] = -sqrt(0.5)*spacing[2]/spacing[0]*direction[0]*direction[2];

  Histogram yHistogramLow(&yBlobs[0], yvec);
  Histogram yHistogramHigh(&yBlobs[1], yvec);
  Histogram *yHistogram[2];
  yHistogram[0] = &yHistogramLow;
  yHistogram[1] = &yHistogramHigh;

  for (int j = 0; j < 2; j++)
    {
    if (!yHistogram[j]->CollapseOne(
          (yPlatePos[j] - origin[1])/spacing[1],
          plateClusterThreshold*xHistogram.GetAverage(), clusterWidthY))
      {
      cerr << "could not find front/back " << j << "\n";;
      foundPlate[2+j] = false;
      continue;
      }
    std::vector<Histogram::Cluster> *yClusters = yHistogram[j]->GetClusters();

    cerr << "PLATE RANGE " << (yPlatePos[j]-origin[1])/spacing[1] << " " << (*yClusters)[0].lowest << " " << (*yClusters)[0].highest << " " << yBlobs[j].size() << "\n";

    for (std::vector<Blob>::iterator it = yBlobs[j].begin();
         it != yBlobs[j].end();
         ++it)
      {
      double y = it->y;
      int yIdx = static_cast<int>(y > 0 ? y + 0.5 : y - 0.5);
      if (yIdx >= (*yClusters)[0].lowest &&
          yIdx <= (*yClusters)[0].highest)
        {
        plateBlobs[2+j].push_back(*it);
        }
      }

    foundPlate[2+j] = plates[2+j].LocateBars(
      &plateBlobs[2+j], xvec, dvec, origin, spacing,
      barSeparation/spacing[0],
      barClusterThreshold, clusterWidth, zLow, zHigh);

    if (foundPlate[2+j])
      {
      cerr << "found end plate " << j << "\n";

      double barVectors[3][3];
      double barCentroids[3][3];
      double barWeights[3];

      for (int k = 0; k < 3; k++)
        {
        FiducialBar *bar = &plates[2+j].GetBars()[k];
        std::vector<Point> *barPoints = bar->GetPoints();
        points->insert(points->end(), barPoints->begin(), barPoints->end());
        barWeights[k] = bar->GetLine(barCentroids[k], barVectors[k]);
        }

      // compute points where diagonal bar intersects vertical bars
      LineIntersection(
        barCentroids[0], barVectors[0], barCentroids[1], barVectors[1],
        corners[4 + 2*j + 0]);
      LineIntersection(
        barCentroids[0], barVectors[0], barCentroids[2], barVectors[2],
        corners[4 + 2*j + 1]);
      }
    }

  // compute the center in data coordinates
  centre[0] = 0.0;
  centre[1] = 0.0;
  centre[2] = 0.0;

  // use side plates
  for (int k = 0; k < 4; k++)
    {
    centre[0] += 0.25*corners[k][0];
    centre[1] += 0.25*corners[k][1];
    centre[2] += 0.25*corners[k][2];
    }

  // improve result with front/back plates if present
  if (foundPlate[2] && foundPlate[3])
    {
    centre[0] *= 0.5;
    centre[1] *= 0.5;
    centre[2] *= 0.5;

    for (int k = 4; k < 8; k++)
      {
      centre[0] += 0.125*corners[k][0];
      centre[1] += 0.125*corners[k][1];
      centre[2] += 0.125*corners[k][2];
      }
    }

  // compute the left-right direction from side plates
  for (int j = 0; j < 3; j++)
    {
    xvec[j] = 0.5*(corners[2][j] + corners[3][j]);
    xvec[j] -= 0.5*(corners[0][j] + corners[1][j]);
    }

  // compute the anterior-posterior direction from front/back plates
  if (foundPlate[2] || foundPlate[3])
    {
    for (int j = 0; j < 3; j++)
      {
      yvec[j] = 0.0;
      if (foundPlate[2])
        {
        yvec[j] += centre[j] - 0.5*(corners[4][j] + corners[5][j]);
        }
      if (foundPlate[3])
        {
        yvec[j] += 0.5*(corners[6][j] + corners[7][j]) - centre[j];
        }
      }
    }
  else
    {
    // if front/back plates are missing, use zvec to compute yvec
    vtkMath::Cross(zvec, xvec, yvec);
    }

  // recompute the z vector from x and y vectors
  vtkMath::Cross(xvec, yvec, zvec);

  // create the frame registration matrix
  BuildMatrix(xvec, yvec, zvec, centre, frameCentre, direction, matrix);

  return true;
}

} // end anonymous namespace

//----------------------------------------------------------------------------
int vtkFrameFinder::FindFrame(
  vtkImageData *image, vtkPolyData *poly,
  const double direction[2], vtkMatrix4x4 *m4x4)
{
  double spacing[3], origin[3];
  int extent[6];
  image->GetSpacing(spacing);
  image->GetOrigin(origin);
  image->GetExtent(extent);

  std::vector<Blob> blobs;
  std::vector<Point> framePoints;

  UpdateBlobs(image, &blobs);

  double matrix[16];
  if (!PositionFrame(&blobs, &framePoints, extent, origin, spacing,
                     direction, matrix))
    {
    cerr << "No frame found!\n";
    }

  m4x4->DeepCopy(matrix);

  if (poly)
    {
    // generate the frame data
    vtkSmartPointer<vtkPoints> points;
    if (poly->GetPoints())
      {
      points = poly->GetPoints();
      }
    else
      {
      points = vtkSmartPointer<vtkPoints>::New();
      }

    vtkSmartPointer<vtkCellArray> cells =
      vtkSmartPointer<vtkCellArray>::New();
    cells->InsertNextCell(0);
    vtkIdType numVerts = 0;

#if 0
    for (std::vector<Blob>::iterator it = blobs.begin();
         it != blobs.end();
         ++it)
      {
      double point[3];
      point[0] = origin[0] + it->x*spacing[0];
      point[1] = origin[1] + it->y*spacing[1];
      point[2] = origin[2] + it->slice*spacing[2];
      vtkIdType ptId = points->InsertNextPoint(point);
      cells->InsertCellPoint(ptId);
      //cout << "blob " << numVerts << " " << point[0] << " " << point[1] << " " << point[2] << "\n";
      numVerts++;
      }

    cells->UpdateCellCount(numVerts);
    poly->SetPoints(points);
    poly->SetVerts(cells);
    }
#else

    for (std::vector<Point>::iterator it = framePoints.begin();
         it != framePoints.end();
         ++it)
      {
      vtkIdType ptId = points->InsertNextPoint(it->x, it->y, it->z);
      cells->InsertCellPoint(ptId);
      numVerts++;
      }

    cells->UpdateCellCount(numVerts);
    poly->SetPoints(points);
    poly->SetVerts(cells);
    }
#endif

  return 1;
}
