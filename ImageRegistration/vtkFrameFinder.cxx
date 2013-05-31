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

#include <math.h>
#include <string.h>

#include <stack>
#include <map>
#include <vector>
#include <utility>

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
vtkImageData *vtkFrameFinder::GetOutput()
{
  return vtkImageData::SafeDownCast(this->GetOutputDataObject(0));
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

//----------------------------------------------------------------------------
namespace {

class lekPointData
{
public:
  lekPointData(int x1, int y1, int z1, int val1) :
    x(x1), y(y1), z(z1), val(val1) {}
  lekPointData() :
    x(0), y(0), z(0), val(0) {}

  int GetX() { return x; }
  int GetY() { return y; }
  int GetZ() { return z; }
  int GetVal() { return val; }

private:
  int x,y,z,val;
};

struct Blob
{
  int numPoints;
  int zSlice;
  double x;
  double y;
  double mass;
};

struct Point
{
  double x;
  double y;
  double z;
};

template <class T>
void GetAndZeroPixel(
  T* array, int x, int y, int z, int dataSizeX, int planeSize, T *value)
{
  *value = array[z*planeSize + y*dataSizeX + x];
  array[z*planeSize + y*dataSizeX + x] = 0;
}

template <class T>
void TestAdd(
  T* array, int x, int y, int z, int dataSizeX, int dataSizeY, int dataSizeZ,
  int threshold, struct Blob *theBlob, std::stack<lekPointData> *basket)
{
  T value;

  if (x >= 0 && y >= 0 && z >= 0 &&
      x < dataSizeX && y < dataSizeY && z < dataSizeZ)
    {
    GetAndZeroPixel(array, x, y, z, dataSizeX, dataSizeY*dataSizeX, &value);
    if (static_cast<int>(value) > threshold)
      {
      lekPointData pD(x, y, z, static_cast<int>(value));
      (*basket).push(pD);
      theBlob->numPoints += 1;
      theBlob->x += value*x;
      theBlob->y += value*y;
      theBlob->mass += value;
      }
    }
}

template <class T>
Blob GetBlob(
  T* array, int x, int y, int z, int dataSizeX, int dataSizeY, int dataSizeZ,
  int threshold)
{
  struct Blob theBlob;

  theBlob.numPoints = 0;
  theBlob.zSlice = z;
  theBlob.x = 0.0;
  theBlob.y = 0.0;
  theBlob.mass = 0.0;
  lekPointData curPt;
  std::stack<lekPointData> basket;

  TestAdd(array, x, y, z, dataSizeX, dataSizeY, dataSizeZ, threshold,
          &theBlob, &basket);

  while (basket.size() > 0 && theBlob.numPoints < 25)
    {
    curPt = basket.top();
    basket.pop();
    TestAdd(array, curPt.GetX() + 1, curPt.GetY(), z,
            dataSizeX, dataSizeY, dataSizeZ,
            threshold, &theBlob, &basket);
    TestAdd(array, curPt.GetX() + 1, curPt.GetY() + 1, z,
            dataSizeX, dataSizeY, dataSizeZ,
            threshold, &theBlob, &basket);
    TestAdd(array, curPt.GetX(), curPt.GetY() + 1, z,
            dataSizeX, dataSizeY, dataSizeZ,
            threshold, &theBlob, &basket);
    TestAdd(array, curPt.GetX() - 1, curPt.GetY() + 1, z,
            dataSizeX, dataSizeY, dataSizeZ,
            threshold, &theBlob, &basket);
    }

  if (theBlob.numPoints > 1 && theBlob.numPoints < 25)
    {
    theBlob.x /= theBlob.mass;
    theBlob.y /= theBlob.mass;
    theBlob.mass /= theBlob.numPoints;
    }
  else
    {
    theBlob.numPoints = 0;
    }

  return theBlob;
}

template <class T>
void ScanVolume(
  std::vector<Blob> *returnList,
  T* inPtr, int dataSizeX, int dataSizeY, int dataSizeZ, int threshold)
{
  int planeSize = dataSizeX*dataSizeY;
  double sum = 0.0;
  double average = 0.0;
  double variance = 0.0;

  T *array = new T [planeSize*dataSizeZ];
  memcpy(array, inPtr, planeSize*dataSizeZ*sizeof(T));

  for (int x = 0; x < dataSizeZ*planeSize; x++)
    {
    sum += array[x];
    }

  average = sum/(planeSize*dataSizeZ);

  for (int x = 0; x < dataSizeZ*planeSize; x++)
    {
    variance += fabs(array[x] - average);
    }
  variance = variance/(planeSize * dataSizeZ);

  threshold = static_cast<int>(average + 2*variance);

  for (int z = 0; z < dataSizeZ; z++)
    {
    for (int y = 0; y < dataSizeY; y++)
      {
      for (int x = 0; x < dataSizeX; x++)
        {
        if (static_cast<int>(array[planeSize*z+dataSizeX*y + x]) > threshold)
          {
          Blob blob = GetBlob(array, x, y, z, dataSizeX, dataSizeY, dataSizeZ,
                              threshold);
          if (blob.numPoints != 0)
            {
            returnList->push_back(blob);
            }
          }
        }
      }
    }

  delete [] array;
}

void FindFrame(
  std::vector<Blob> *returnList, void *address, const int extent[6],
  int scalarType, int threshold)
{
  int dataSizeX = extent[1] - extent[0] + 1;
  int dataSizeY = extent[3] - extent[2] + 1;
  int dataSizeZ = extent[5] - extent[4] + 1;

  switch (scalarType)
    {
    case VTK_DOUBLE:
      ScanVolume(returnList, static_cast<double *>(address),
                 dataSizeX, dataSizeY, dataSizeZ, threshold);
      break;
    case VTK_FLOAT:
      ScanVolume(returnList, static_cast<float *>(address),
                 dataSizeX, dataSizeY, dataSizeZ, threshold);
      break;
    case VTK_LONG:
      ScanVolume(returnList, static_cast<long *>(address),
                 dataSizeX, dataSizeY, dataSizeZ, threshold);
      break;
    case VTK_UNSIGNED_LONG:
      ScanVolume(returnList, static_cast<unsigned long *>(address),
                 dataSizeX, dataSizeY, dataSizeZ, threshold);
      break;
    case VTK_INT:
      ScanVolume(returnList, static_cast<int *>(address),
                 dataSizeX, dataSizeY, dataSizeZ, threshold);
      break;
    case VTK_UNSIGNED_INT:
      ScanVolume(returnList, static_cast<unsigned int *>(address),
                 dataSizeX, dataSizeY, dataSizeZ, threshold);
      break;
    case VTK_SHORT:
      ScanVolume(returnList, static_cast<short *>(address),
                 dataSizeX, dataSizeY, dataSizeZ, threshold);
      break;
    case VTK_UNSIGNED_SHORT:
      ScanVolume(returnList, static_cast<unsigned short *>(address),
                 dataSizeX, dataSizeY, dataSizeZ, threshold);
      break;
    case VTK_CHAR:
      ScanVolume(returnList, static_cast<char *>(address),
                 dataSizeX, dataSizeY, dataSizeZ, threshold);
      break;
    case VTK_UNSIGNED_CHAR:
      ScanVolume(returnList, static_cast<unsigned char *>(address),
                 dataSizeX, dataSizeY, dataSizeZ, threshold);
      break;
    case VTK_SIGNED_CHAR:
      ScanVolume(returnList, static_cast<signed char *>(address),
                 dataSizeX, dataSizeY, dataSizeZ, threshold);
      break;
    }
}

void UpdateBlobs(std::vector<Blob> *blobs, vtkImageData *input)
{
  FindFrame(blobs, input->GetScalarPointer(), input->GetExtent(),
            input->GetScalarType(), 15);
}

class Histogram
{
public:
  Histogram(std::map<int, double> *dict);
  double ReturnCluster(int *leftBinP, int *rightBinP);
  bool Collapse(double distance);

  class Bin
  {
  public:
    Bin(int l, int r, double s) :
      leftBin(l), rightBin(r), sum(s) {}

    int leftBin;
    int rightBin;
    double sum;
  };

  std::vector<Bin> *GetBins() { return &this->Bins; }
  std::map<int, double> *GetRawBins() { return this->HistoBins; }

private:
  std::vector<Bin> Bins;
  std::map<int, double> *HistoBins;
  int MaxPeakWidth;
  double average;
  double biggest;
  double smallest;
  double KeepThreshold;
};

Histogram::Histogram(std::map<int, double> *dict)
{
  double sum = 0.0;
  double maxval = 0.0;
  double minval = 1e30;

  for (std::map<int, double>::iterator it = dict->begin();
       it != dict->end();
       ++it)
    {
    double v = it->second;
    sum += v;
    maxval = (v < maxval ? maxval : v);
    minval = (v > minval ? minval : v);
    }

  this->MaxPeakWidth = 20; // 10; must be divided by spacing
  this->HistoBins = dict;
  this->average = sum/dict->size();
  this->biggest = maxval;
  this->smallest = minval;
  this->KeepThreshold = 0.5*this->average/4; // multiply by x*y spacing
}

double Histogram::ReturnCluster(int *leftBinP, int *rightBinP)
{
  // Picks biggest peak in histogram, returns tuple consisting
  // of range of peak and area of peak.
  // Removes bins from HistoBins
  // Multiple calls return the multiple biggest peaks in the
  // histogram

  // first find the peak bin
  std::map<int, double> *dict = this->HistoBins;
  int maxBin = 0;
  double maxBinValue = 0;
  for (std::map<int, double>::iterator it = dict->begin();
       it != dict->end();
       ++it)
    {
    int i = it->first;
    double v = it->second;
    if (v > maxBinValue)
      {
      maxBin = i;
      maxBinValue= v;
      }
    }
  double sum = maxBinValue;
  dict->erase(maxBin);

  // then count down from the middle
  int leftBin = maxBin;
  for (;;)
    {
    std::map<int, double>::iterator it = dict->find(leftBin - 1);
    if (it == dict->end()) { break; }
    double v = it->second;
    if (v <= this->KeepThreshold) { break; }
    sum += v;
    dict->erase(it);
    leftBin--;
    }

  // then count up from the middle
  int rightBin = maxBin;
  for (;;)
    {
    std::map<int, double>::iterator it = dict->find(rightBin + 1);
    if (it == dict->end()) { break; }
    double v = it->second;
    if (v <= this->KeepThreshold) { break; }
    sum += v;
    dict->erase(it);
    rightBin++;
    }

  *leftBinP = leftBin;
  *rightBinP = rightBin;
  return sum;
}

bool Histogram::Collapse(double distance)
{
  // distance = ideal distance between histogram peaks in
  // data coordinates
  std::vector<Bin> *bins = &this->Bins;
  std::map<int, double> *dict = this->HistoBins;

  bins->clear();
  while (!dict->empty())
    {
    int leftBin, rightBin;
    double sum = this->ReturnCluster(&leftBin, &rightBin);
    if (rightBin - leftBin <= this->MaxPeakWidth &&
        sum > this->average)
      {
      cerr << rightBin << " - " << leftBin << " = " << this->MaxPeakWidth << " : " << sum << " > " << this->average << "\n";
      bins->push_back(Bin(leftBin, rightBin, sum));
      }
    }

  int bestDist = 100;
  std::vector<Bin>::iterator bestBin1 = bins->end();
  std::vector<Bin>::iterator bestBin2 = bins->end();

  for (std::vector<Bin>::iterator it = bins->begin();
       it != bins->end();
       ++it)
    {
    int binPos = (it->leftBin + it->rightBin)/2;
    for (std::vector<Bin>::iterator jt = bins->begin();
         jt != bins->end();
         ++jt)
      {
      int binPos2 = (jt->leftBin + jt->rightBin)/2;
      int curDist = binPos - binPos2;
      curDist = (curDist >= 0 ? curDist : -curDist);
      curDist -= distance;
      curDist = (curDist >= 0 ? curDist : -curDist);
      if (curDist < bestDist)
        {
        bestDist = curDist;
        bestBin1 = it;
        bestBin2 = jt;
        }
      }
    }

  if (bestBin1 == bins->end() || bestBin2 == bins->end())
    {
    return false;
    }

  if (bestBin1->leftBin > bestBin2->leftBin)
    {
    std::vector<Bin>::iterator tmp = bestBin2;
    bestBin2 = bestBin1;
    bestBin1 = tmp;
    }

  cerr << "bestDist " << bestDist << "\n";
  cerr << "pos " << (bestBin1->leftBin + bestBin1->rightBin)/2
       << " " << (bestBin2->leftBin + bestBin2->rightBin)/2
       << " sum " << bestBin1->sum << " " << bestBin2->sum << "\n";

  Bin bin1 = *bestBin1;
  Bin bin2 = *bestBin2;

  bins->clear();
  bins->push_back(bin1);
  bins->push_back(bin2);

  return true;
}

// collapse the blob points into bins by x or y value
void MakeHistogram(
  std::map<int, double> *dict, std::vector<Blob> *blobs, int i)
{
  for (std::vector<Blob>::iterator it = blobs->begin();
       it != blobs->end();
       ++it)
    {
    int j = static_cast<int>(i == 0 ? it->x : it->y);
    std::map<int, double>::iterator jt = dict->find(j);
    if (jt == dict->end())
      {
      dict->insert(std::make_pair(j, it->mass));
      }
    else
      {
      jt->second += it->mass;
      }
    }
}

class FiducialBar
{
public:
  FiducialBar() {
    this->points = 0;
    this->eigenvector[0] = this->eigenvector[1] = this->eigenvector[2] = 0;
    this->centreOfMass[0] = this->centreOfMass[1] = this->centreOfMass[2] = 0;
    this->eigenvalue = 0; }

  void SetPoints(std::vector<Point> *p) {
    this->points = p;
    this->centreOfMass[0] = 0;
    this->centreOfMass[1] = 0;
    this->centreOfMass[2] = 0; }

  std::vector<Point> *GetPoints() {
    return this->points; }

  void ClipFiducialEnds(double zMin, double zMax);

  void ExtractLinesFromPoints();

  double CullOutliersAndReturnRMS(double maxDist);

  void GetLine(double p[3], double v[3]);

private:
  void ComputeCentreOfMass(double com[3]);

  std::vector<Point> *points;
  double eigenvalue;
  double eigenvector[3];
  double centreOfMass[3];
};

void FiducialBar::GetLine(double p[3], double v[3])
{
  p[0] = this->centreOfMass[0];
  p[1] = this->centreOfMass[1];
  p[2] = this->centreOfMass[2];

  v[0] = this->eigenvector[0];
  v[1] = this->eigenvector[1];
  v[2] = this->eigenvector[2];

  if (v[2] < 0)
    {
    v[0] = -v[0];
    v[1] = -v[1];
    v[2] = -v[2];
    }
}

void FiducialBar::ClipFiducialEnds(double minZ, double maxZ)
{
  size_t j = 0;
  for (size_t i = 0; i < this->points->size(); i++)
    {
    Point &pi = (*this->points)[i];
    if (pi.z > minZ && pi.z < maxZ)
      {
      (*this->points)[j++] = pi;
      }
    }
  this->points->resize(j);
}

void FiducialBar::ComputeCentreOfMass(double com[3])
{
  double count = 0;
  com[0] = 0.0;
  com[1] = 0.0;
  com[2] = 0.0;

  for (std::vector<Point>::iterator it = this->points->begin();
       it < this->points->end();
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

void FiducialBar::ExtractLinesFromPoints()
{
  double C[3][3] = { {0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0} };
  this->ComputeCentreOfMass(this->centreOfMass);
  double *Q = this->centreOfMass;
  for (std::vector<Point>::iterator it = this->points->begin();
       it < this->points->end();
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

  this->eigenvalue = eigenvalues[0];
  this->eigenvector[0] = eigenvectors[0][0];
  this->eigenvector[1] = eigenvectors[1][0];
  this->eigenvector[2] = eigenvectors[2][0];

  double *v = this->eigenvector;
  cerr << "(" << Q[0] << ", " << Q[1] << ", " << Q[2] << ") ";
  cerr << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")\n";
}

double FiducialBar::CullOutliersAndReturnRMS(double maxDist)
{
  // Check each point in the bar for distance from the extracted line
  // If the distance exceeds maxDist, remove the point from the list.

  double maxDistSquared = maxDist*maxDist;
  double sum = 0.0;
  double *com = this->centreOfMass;

  size_t j = 0;
  for (size_t i = 0; i < this->points->size(); i++)
    {
    Point &pi = (*this->points)[i];

    // compute distance from point to line with pythagoras
    double p[3];
    double *v = this->eigenvector;
    p[0] = pi.x - com[0];
    p[1] = pi.y - com[1];
    p[2] = pi.z - com[2];
    double dotpp = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
    double dotpv = v[0]*p[0] + v[1]*p[1] + v[2]*p[2];
    double distSquared = dotpp - dotpv*dotpv;

    if (distSquared < maxDistSquared)
      {
      sum += distSquared;
      (*this->points)[j++] = pi;
      }
    }

  // points too far from line have been discarded
  this->points->resize(j);

  // compute the RMS
  if (j != 0)
    {
    sum = sqrt(sum/static_cast<double>(j));
    }

  return sum;
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

bool PositionFrame(
  std::vector<Blob> *blobs, std::vector<Point> *points,
  const double origin[3], const double spacing[3],
  const double direction[3], double matrix[16])
{
  std::map<int, double> xDict;

  cerr << "spacing " << spacing[0] << " " << spacing[1] << " " << spacing[2] << "\n";

  MakeHistogram(&xDict, blobs, 0);
  Histogram xHistogram(&xDict);
  if (!xHistogram.Collapse(190.0/spacing[0]))
    {
    cerr << "could not find sides\n";
    return false;
    }
  std::vector<Histogram::Bin> *xBins = xHistogram.GetBins();

  int zMin = 10000;
  int zMax = -1;

  std::vector<Blob> lowX;
  std::vector<Blob> highX;

  for (std::vector<Blob>::iterator it = blobs->begin();
       it != blobs->end();
       ++it)
    {
    zMin = (zMin < it->zSlice ? zMin : it->zSlice);
    zMax = (zMax > it->zSlice ? zMax : it->zSlice);
    if (static_cast<int>(it->x) >= (*xBins)[0].leftBin &&
        static_cast<int>(it->x) <= (*xBins)[0].rightBin)
      {
      lowX.push_back(*it);
      }
    else if (static_cast<int>(it->x) >= (*xBins)[1].leftBin &&
             static_cast<int>(it->x) <= (*xBins)[1].rightBin)
      {
      highX.push_back(*it);
      }
    }

  // the range of z values to actually use
  int zLow = zMin + (zMax - zMin)/4;
  int zHigh = zMax - (zMax - zMin)/4;

  std::map<int, double> lowBlobs;
  std::map<int, double> highBlobs;

  MakeHistogram(&lowBlobs, &lowX, 1);
  Histogram lowXHisto(&lowBlobs);
  if (!lowXHisto.Collapse(120.0/spacing[1]))
    {
    cerr << "could not find left\n";
    return false;
    }
  std::vector<Histogram::Bin> *lowXBins = lowXHisto.GetBins();

  MakeHistogram(&highBlobs, &highX, 1);
  Histogram highXHisto(&highBlobs);
  if (!highXHisto.Collapse(120.0/spacing[1]))
    {
    cerr << "could not find right\n";
    return false;
    }
  std::vector<Histogram::Bin> *highXBins = highXHisto.GetBins();

  double corners[4][3];

  for (int j = 0; j < 2; j++)
  {

  std::vector<Histogram::Bin> *barBins =
    (j == 0 ? lowXBins : highXBins);
  std::vector<Blob> *barBlobs =
    (j == 0 ? &lowX : &highX);

  std::vector<Point> firstBarPoints;
  std::vector<Point> lastBarPoints;
  std::vector<Point> diagonalBarPoints;

  for (std::vector<Blob>::iterator it = barBlobs->begin();
       it != barBlobs->end();
       ++it)
    {
    if (it->zSlice < zLow || it->zSlice > zHigh)
      {
      continue;
      }

    int idx = static_cast<int>(it->y);

    Point p;
    p.x = origin[0] + spacing[0]*it->x;
    p.y = origin[1] + spacing[1]*it->y;
    p.z = origin[2] + spacing[2]*it->zSlice;

    if (idx >= (*barBins)[0].leftBin && idx <= (*barBins)[0].rightBin)
      {
      firstBarPoints.push_back(p);
      points->push_back(p);
      }
    else if (idx >= (*barBins)[1].leftBin && idx <= (*barBins)[1].rightBin)
      {
      lastBarPoints.push_back(p);
      points->push_back(p);
      }
    else if (idx > (*barBins)[0].rightBin && idx < (*barBins)[1].leftBin)
      {
      diagonalBarPoints.push_back(p);
      points->push_back(p);
      }
    }

  FiducialBar firstBar;
  firstBar.SetPoints(&firstBarPoints);
  firstBar.ExtractLinesFromPoints();
  firstBar.CullOutliersAndReturnRMS(15);
  firstBar.ExtractLinesFromPoints();
  double firstCentre[3], firstVector[3];
  firstBar.GetLine(firstCentre, firstVector);

  FiducialBar lastBar;
  lastBar.SetPoints(&lastBarPoints);
  lastBar.ExtractLinesFromPoints();
  lastBar.CullOutliersAndReturnRMS(15);
  lastBar.ExtractLinesFromPoints();
  double lastCentre[3], lastVector[3];
  lastBar.GetLine(lastCentre, lastVector);

  FiducialBar diagonalBar;
  diagonalBar.SetPoints(&diagonalBarPoints);
  diagonalBar.ExtractLinesFromPoints();
  diagonalBar.CullOutliersAndReturnRMS(15);
  diagonalBar.ExtractLinesFromPoints();
  double diagonalCentre[3], diagonalVector[3];
  diagonalBar.GetLine(diagonalCentre, diagonalVector);

  double *p1 = corners[2*j+0];
  double *p2 = corners[2*j+1];
  LineIntersection(
    firstCentre, firstVector, diagonalCentre, diagonalVector, p1);
  LineIntersection(
    lastCentre, lastVector, diagonalCentre, diagonalVector, p2);

  cerr << "p1 = " << p1[0] << " " << p1[1] << " " << p1[2] << "\n";
  cerr << "p2 = " << p2[0] << " " << p2[1] << " " << p2[2] << "\n";
  }

  // compute the horizontal direction
  double xvec[3];
  for (int j = 0; j < 3; j++)
    {
    xvec[j] = 0.5*(corners[2][j] + corners[3][j]);
    xvec[j] -= 0.5*(corners[0][j] + corners[1][j]);
    }
  vtkMath::Normalize(xvec);

  // compute the diagonal direction
  double dvec[3];
  for (int j = 0; j < 3; j++)
    {
    dvec[j] = 0.5*(corners[1][j] + corners[3][j]);
    dvec[j] -= 0.5*(corners[0][j] + corners[2][j]);
    }
  vtkMath::Normalize(dvec);

  // orthogonalize wrt horizontal
  double tvec[3];
  vtkMath::Cross(dvec, xvec, tvec);
  vtkMath::Cross(xvec, tvec, dvec);
  vtkMath::Normalize(dvec);

  // get the yvec and zvec from dvec
  double yvec[3];
  yvec[0] = dvec[0] + tvec[0];
  yvec[1] = dvec[1] + tvec[1];
  yvec[2] = dvec[2] + tvec[2];
  vtkMath::Normalize(yvec);

  double zvec[3];
  zvec[0] = dvec[0] - tvec[0];
  zvec[1] = dvec[1] - tvec[1];
  zvec[2] = dvec[2] - tvec[2];
  vtkMath::Normalize(zvec);

  // compute the center in data coordinates
  double centre[3];
  centre[0] = 0.0;
  centre[1] = 0.0;
  centre[2] = 0.0;
  for (int k = 0; k < 4; k++)
    {
    centre[0] += 0.25*corners[k][0];
    centre[1] += 0.25*corners[k][1];
    centre[2] += 0.25*corners[k][2];
    }

  cerr << "centre = " << centre[0] << " " << centre[1] << " " << centre[2] << "\n";

  double dsign[3];
  dsign[0] = (direction[0] < 0 ? -1.0 : +1.0);
  dsign[1] = (direction[1] < 0 ? -1.0 : +1.0);
  dsign[2] = (direction[2] < 0 ? -1.0 : +1.0);

  // build the matrix
  for (int k = 0; k < 3; k++)
    {
    matrix[4*k + 0] = xvec[k]*dsign[0];
    matrix[4*k + 1] = yvec[k]*dsign[1];
    matrix[4*k + 2] = zvec[k]*dsign[2];
    matrix[12 + k] = 0.0;
    }

  double tc[3];
  tc[0] = matrix[0]*centre[0] + matrix[1]*centre[1] + matrix[2]*centre[2];
  tc[1] = matrix[4]*centre[0] + matrix[5]*centre[1] + matrix[6]*centre[2];
  tc[2] = matrix[8]*centre[0] + matrix[9]*centre[1] + matrix[10]*centre[2];

  matrix[3] = 100.0 - tc[0];
  matrix[7] = 100.0 - tc[1];
  matrix[11] = 100.0 - tc[2];
  matrix[15] = 1.0;

  return true;
}

} // end anonymous namespace

//----------------------------------------------------------------------------
int vtkFrameFinder::FindFrame(
  vtkImageData *image, vtkPolyData *poly,
  const double direction[2], vtkMatrix4x4 *m4x4)
{
  double spacing[3], origin[3];
  image->GetSpacing(spacing);
  image->GetOrigin(origin);

  std::vector<Blob> blobs;
  std::vector<Point> framePoints;

  UpdateBlobs(&blobs, image);

  double matrix[16];
  if (!PositionFrame(&blobs, &framePoints,
      origin, spacing, direction, matrix))
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

    /*
    for (std::vector<Blob>::iterator it = blobs.begin();
         it != blobs.end();
         ++it)
      {
      double point[3];
      point[0] = origin[0] + it->x*spacing[0];
      point[1] = origin[1] + it->y*spacing[1];
      point[2] = origin[2] + it->zSlice*spacing[2];
      vtkIdType ptId = points->InsertNextPoint(point);
      cells->InsertCellPoint(ptId);
      //cout << "blob " << numVerts << " " << point[0] << " " << point[1] << " " << point[2] << "\n";
      numVerts++;
      }

    cells->UpdateCellCount(numVerts);
    poly->SetPoints(points);
    poly->SetVerts(cells);
    }
    */

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

  return 1;
}
