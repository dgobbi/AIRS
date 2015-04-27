/*=========================================================================

Copyright (c) 2004 Atamai, Inc.

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
#include "vtkImageMRIBrainExtractor.h"
#include "vtkObjectFactory.h"
#include "vtkImageData.h"
#include "vtkPolyDataToImageStencil.h"
#include "vtkImageStencil.h"
#include "vtkPolyData.h"
#include "vtkSphereSource.h"
#include "vtkLinearSubdivisionFilter.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkIdList.h"
#include "vtkIdListCollection.h"
#include "vtkPolyDataNormals.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkMath.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkImageIterator.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkVersion.h"

#include <math.h>
#include <iostream>

#include <vector>
#include <algorithm>
#include <numeric>

// A macro to assist VTK 5 backwards compatibility
#if VTK_MAJOR_VERSION >= 6
#define SET_INPUT_DATA SetInputData
#define SET_SOURCE_DATA SetSourceData
#define SET_STENCIL_DATA SetStencilData
#else
#define SET_INPUT_DATA SetInput
#define SET_SOURCE_DATA SetSource
#define SET_STENCIL_DATA SetStencil
#endif

namespace {

// A convenience class for using C arrays in STL vectors
class myVertIds
{
public:
  int abc[3];

  myVertIds() {};
  myVertIds(const int vIds[3]) {
    abc[0] = vIds[0]; abc[1] = vIds[1]; abc[2] = vIds[2]; };
  myVertIds(const myVertIds &vIds) {
    abc[0] = vIds[0]; abc[1] = vIds[1]; abc[2] = vIds[2]; };

  const int& operator[](int i) const { return abc[i]; };
  int& operator[](int i) { return abc[i]; };
};

// A class to avoid the overhead of small STL vectors
class myPoint
{
public:
  double xyz[3];

  myPoint() {};
  myPoint(const double point[3]) {
    xyz[0] = point[0]; xyz[1] = point[1]; xyz[2] = point[2]; };
  myPoint(const float point[3]) {
    xyz[0] = point[0]; xyz[1] = point[1]; xyz[2] = point[2]; };
  myPoint(const myPoint &point) {
    xyz[0] = point[0]; xyz[1] = point[1]; xyz[2] = point[2]; };

  const double& operator[](int i) const { return xyz[i]; };
  double& operator[](int i) { return xyz[i]; };
};

}

vtkStandardNewMacro(vtkImageMRIBrainExtractor);

//----------------------------------------------------------------------------
vtkImageMRIBrainExtractor::vtkImageMRIBrainExtractor()
{
  this->BrainMesh = vtkPolyData::New();
  // Defaults
  this->BT = 0.5;
  this->NumberOfIterations = 1000;
  this->NumberOfTessellations = 4;
  this->D1 = 20.0; //mm
  this->D2 = 10.0; //mm
  this->RMin = 3.0; //mm
  this->RMax = 10.0; //mm

  this->BrainExtent[0] = VTK_INT_MIN;
  this->BrainExtent[1] = VTK_INT_MAX;
  this->BrainExtent[2] = VTK_INT_MIN;
  this->BrainExtent[3] = VTK_INT_MAX;
  this->BrainExtent[4] = VTK_INT_MIN;
  this->BrainExtent[5] = VTK_INT_MAX;
}

//----------------------------------------------------------------------------
vtkImageMRIBrainExtractor::~vtkImageMRIBrainExtractor()
{
  this->BrainMesh->Delete();
}

//----------------------------------------------------------------------------
vtkPolyData *vtkImageMRIBrainExtractor::GetBrainMesh()
{
  return this->BrainMesh;
}

//----------------------------------------------------------------------------
void vtkImageMRIBrainExtractor::ComputeMeshCentroid(
  vtkPolyData *data, double cen[3])
{
  vtkPoints *points = data->GetPoints();
  vtkCellArray *polys = data->GetPolys();
  double vol = 0.0;
  cen[0] = 0.0;
  cen[1] = 0.0;
  cen[2] = 0.0;

  /* only use polys, ignore strips */
  if (polys)
    {
    vtkIdType n = polys->GetNumberOfCells();
    vtkIdType l = 0;
    vtkIdType *ptIds;
    vtkIdType nPts;
    for (vtkIdType i = 0; i < n; i++)
      {
      polys->GetCell(l, nPts, ptIds);
      l += nPts + 1;
      double v1[3];
      double v2[3];
      double p0[3];
      double p1[3];
      vtkIdType m = nPts - 1;
      points->GetPoint(ptIds[m], p1);
      points->GetPoint(ptIds[0], p0);
      v2[0] = p1[0] - p0[0];
      v2[1] = p1[1] - p0[1];
      v2[2] = p1[2] - p0[2];

      for (vtkIdType j = 1; j < m; j++)
        {
        v1[0] = v2[0];
        v1[1] = v2[1];
        v1[2] = v2[2];

        points->GetPoint(ptIds[j], p1);
        v2[0] = p1[0] - p0[0];
        v2[1] = p1[1] - p0[1];
        v2[2] = p1[2] - p0[2];

        /* scalar triple product gives volume times six*/
        double pvol = 0.0;
        pvol += p0[0]*(v1[1] * v2[2] - v1[2] * v2[1]);
        pvol += p0[1]*(v1[2] * v2[0] - v1[0] * v2[2]);
        pvol += p0[2]*(v1[0] * v2[1] - v1[1] * v2[0]);
        vol += pvol;

        /* this gives tetrahedron centroid times 4*6 */
        cen[0] += pvol*(3*p0[0] + v1[0] + v2[0]);
        cen[1] += pvol*(3*p0[1] + v1[1] + v2[1]);
        cen[2] += pvol*(3*p0[2] + v1[2] + v2[2]);
        }
      }
    }

  if (vol != 0)
    {
    cen[0] /= 4*vol;
    cen[1] /= 4*vol;
    cen[2] /= 4*vol;
    }
}

namespace {
//----------------------------------------------------------------------------
// Description:
// This templated function executes the filter for any type of data.
template <class IT>
static void vtkBECalculateInitialParameters(
  vtkImageData *inData, IT *vtkNotUsed(inPtr), int extent[6],
  double &T2, double &T98, double &TH,
  double &Tm, double COG[3], double &R)
{
  int scalarTypeMin = static_cast<int>(inData->GetScalarTypeMin());
  int scalarTypeMax = static_cast<int>(inData->GetScalarTypeMax());

  // Initialize the values to output
  T2 = scalarTypeMin;
  Tm = scalarTypeMin + 0.5;
  T98 = scalarTypeMin + 1.0;
  COG[0] = COG[1] = COG[2] = 0.0;
  R = 1.0;

  int nBins = scalarTypeMax - scalarTypeMin + 1;
  size_t *hist = new size_t[nBins];

  // initialize the histogram
  for (int j = 0; j < nBins; j++)
    {
    hist[j] = 0;
    }

  size_t voxelCount = (extent[1] - extent[0] + 1);
  voxelCount *= (extent[3] - extent[2] + 1);
  voxelCount *= (extent[5] - extent[4] + 1);

  size_t lowerThreshold = static_cast<size_t>(0.02*voxelCount + 0.5);
  size_t upperThreshold = static_cast<size_t>(0.98*voxelCount + 0.5);

  // accumulate histogram bins
  vtkImageIterator<IT> inIter(inData, extent);
  while (!inIter.IsAtEnd())
    {
    IT *tmpPtr = inIter.BeginSpan();
    IT *tmpPtrEnd = inIter.EndSpan();
    while (tmpPtr != tmpPtrEnd)
      {
      int idx = static_cast<int>(*tmpPtr - scalarTypeMin);
      hist[idx]++;
      tmpPtr++;
      }
    inIter.NextSpan();
    }

  // compute thresholds
  size_t histogramSum = 0;

  for (int bin = 0; bin < nBins; bin++)
    {
    size_t f = hist[bin];
    histogramSum += f;
    int v = bin + scalarTypeMin;

    if (histogramSum <= lowerThreshold)
      {
      T2 = v;
      }
    else if (histogramSum <= upperThreshold)
      {
      T98 = v;
      }
    }

  TH = T2 + (0.10*(T98-T2));

  //------------------------
  double XMoment = 0.0;
  double YMoment = 0.0;
  double ZMoment = 0.0;
  double totalMass = 0.0;

  double *spacing = inData->GetSpacing();
  double *origin = inData->GetOrigin();

  vtkIdType count = 0;
  vtkIdType inIncX, inIncY, inIncZ;
  inData->GetContinuousIncrements(extent, inIncX, inIncY, inIncZ);
  IT *tmpPtr = static_cast<IT *>(inData->GetScalarPointerForExtent(extent));

  for (int idx2 = extent[4]; idx2 <= extent[5]; idx2++)
    {
    for (int idx1 = extent[2]; idx1 <= extent[3]; idx1++)
      {
      for (int idx0 = extent[0]; idx0 <= extent[1]; idx0++)
        {
        double mass = static_cast<double>(*tmpPtr);
        if (mass > TH)
          {
          // Limit our mass so it's not an outlier value
          mass = ((mass <= T98) ? mass : T98);
          XMoment += mass*idx0;
          YMoment += mass*idx1;
          ZMoment += mass*idx2;
          totalMass += mass;
          count += 1;
          }
        tmpPtr++;
        }
      tmpPtr += inIncY;
      }
    tmpPtr += inIncZ;
    }

  if (totalMass == 0)
    {
    vtkGenericWarningMacro("In vtkMRIBrainExtractor, image is all black");
    return;
    }

  COG[0] = ((XMoment / totalMass) * spacing[0]) + origin[0];
  COG[1] = ((YMoment / totalMass) * spacing[1]) + origin[1];
  COG[2] = ((ZMoment / totalMass) * spacing[2]) + origin[2];

  double voxelVolume = fabs(spacing[0]*spacing[1]*spacing[2]);
  double totalVolume = voxelVolume*count;

  R = pow( (3.0*totalVolume)/(4.0*vtkMath::Pi()), 1.0/3.0 );

  // Compute the median value within a sphere of radius R
  tmpPtr = static_cast<IT *>(inData->GetScalarPointerForExtent(extent));

  // initialize the histogram
  for (int j = 0; j < nBins; j++)
    {
    hist[j] = 0;
    }
  voxelCount = 0;

  double R2 = R*R;
  for (int idx2 = extent[4]; idx2 <= extent[5]; idx2++)
    {
    double zz = idx2*spacing[2] + origin[2] - COG[2];
    zz *= zz;
    for (int idx1 = extent[2]; idx1 <= extent[3]; idx1++)
      {
      double yy = idx1*spacing[1] + origin[1] - COG[1];
      yy = zz + yy*yy;
      double xx = extent[0]*spacing[0] + origin[0] - COG[0];
      for (int idx0 = extent[0]; idx0 <= extent[1]; idx0++)
        {
        if (T2 < *tmpPtr && *tmpPtr < T98)
          {
          double distance2 = yy + xx*xx;

          // Are we within R?
          if (distance2 < R2)
            {
            hist[static_cast<int>(*tmpPtr - scalarTypeMin)] += 1;
            voxelCount++;
            }
          }
        tmpPtr++;
        xx += spacing[0];
        }
      tmpPtr += inIncY;
      }
    tmpPtr += inIncZ;
    }

  histogramSum = 0;
  size_t medianThreshold = voxelCount/2;
  for (int bin = 0; bin < nBins; bin++)
    {
    int f = hist[bin];
    histogramSum += f;

    if (histogramSum > medianThreshold)
      {
      Tm = static_cast<double>(bin) -
        (histogramSum - 0.5*voxelCount)/f +
        scalarTypeMin;
      break;
      }
    }

  delete [] hist;
}

//----------------------------------------------------------------------------
// Description:
static void vtkBEBuildAndLinkPolyData(
  double COG[3], double R, int Nsubs,
  vtkPolyData *brainPolyData,
  vtkIdListCollection *pointNeighbourList)
{
  int nPoints, ptId, cellIdx, cellId, ptIdx2, ptId2;
  vtkIdType *pointCells, npts, *pts;
  unsigned short nCells;

  // Icosahedron - a 20-sided polygon
  vtkSphereSource *icosahedron = vtkSphereSource::New();
  icosahedron->SetCenter(COG[0], COG[1], COG[2]);
  icosahedron->SetRadius(R);
  icosahedron->SetPhiResolution(4);
  icosahedron->SetThetaResolution(5);
  icosahedron->Update();

  // Subdivide each triangle into 4.
  vtkLinearSubdivisionFilter *subdivideSphere =
    vtkLinearSubdivisionFilter::New();
  subdivideSphere->SET_INPUT_DATA(icosahedron->GetOutput());
  subdivideSphere->SetNumberOfSubdivisions(Nsubs);
  subdivideSphere->Update();

  // Contraint sphere for smoothing
  vtkSphereSource *constraintSphere = vtkSphereSource::New();
  constraintSphere->SetCenter(COG[0], COG[1], COG[2]);
  constraintSphere->SetRadius(R);
  constraintSphere->SetPhiResolution(100);
  constraintSphere->SetThetaResolution(100);
  constraintSphere->Update();

  // Smooth the subdivided sphere
  vtkSmoothPolyDataFilter *smoothSphere = vtkSmoothPolyDataFilter::New();
  smoothSphere->SET_INPUT_DATA(subdivideSphere->GetOutput());
  smoothSphere->SET_SOURCE_DATA(constraintSphere->GetOutput());
  smoothSphere->Update();

  // The brain sphere
  brainPolyData->DeepCopy(smoothSphere->GetOutput());
  brainPolyData->BuildLinks();

  nPoints = brainPolyData->GetNumberOfPoints();

  // Create a list of neighbours for each point
  for (ptId = 0; ptId < nPoints ; ptId++)
    {
    // Create a new vtkIdList to hold the Ids for ptId's neighbours
    vtkIdList *myNeighbours = vtkIdList::New();

    // The cells attached to the target point
    brainPolyData->GetPointCells( ptId, nCells, pointCells );

    // The points inside the attached cells
    for (cellIdx = 0; cellIdx < nCells; cellIdx++)
      {
      cellId = pointCells[cellIdx];
      brainPolyData->GetCellPoints(cellId, npts, pts);

      // If the point is not our target point, it's a neighbour
      for ( ptIdx2 = 0; ptIdx2 < npts ; ptIdx2++ )
        {
        ptId2 = pts[ptIdx2];
        if (ptId != ptId2)
          {
          myNeighbours->InsertUniqueId(ptId2);
          }
        }
      }

    pointNeighbourList->AddItem( myNeighbours );
    myNeighbours->Delete();
    }

  // Clean up
  icosahedron->Delete();
  subdivideSphere->Delete();
  constraintSphere->Delete();
  smoothSphere->Delete();
}

//----------------------------------------------------------------------------
// This templated function executes the filter for any type of data.
template <class IT>
void vtkImageMRIBrainExtractorExecute(
  vtkImageMRIBrainExtractor *self,
  vtkImageData *inData, IT *inPtr)
{
  // The parameters for the algorithm
  double T2, T98, TH, Tm;
  double COG[3], R;

  int extent[6];
  double origin[3], spacing[3];
  inData->GetExtent(extent);
  inData->GetSpacing(spacing);
  inData->GetOrigin(origin);

  int brainExtent[6];
  self->GetBrainExtent(brainExtent);
  for (int k = 0; k < 3; k++)
    {
    if (brainExtent[2*k] < extent[2*k])
      {
      brainExtent[2*k] = extent[2*k];
      }
    if (brainExtent[2*k+1] > extent[2*k+1])
      {
      brainExtent[2*k+1] = extent[2*k+1];
      }
    }

  // vtkImageData-based parameters
  vtkBECalculateInitialParameters(
    inData, inPtr, brainExtent, T2, T98, TH, Tm, COG, R);

  // vtkPolyData time
  vtkPolyData *brainPolyData = vtkPolyData::New();
  vtkIdListCollection *pointNeighbourList = vtkIdListCollection::New();

  // Initialize a sphere inside the brain and for each point on the sphere,
  // build a list of its first-order neighbours
  double RSphere = 0.5*R;

  vtkBEBuildAndLinkPolyData(
    COG, RSphere, self->GetNumberOfTessellations(),
    brainPolyData, pointNeighbourList);

  vtkIdType nPoints = brainPolyData->GetNumberOfPoints();

  int nIterations, iteration, nNeighbours;
  double dotp, suml2, this_suml2;
  std::vector<double> Id1;
  std::vector<double> Id2;
  double value, d;
  double direction[3];
  int end[3], start[3];
  double fend[3], fstart[3];
  myPoint newtarget;
  double location[3];

  // BET variables
  double n_hat[3], s[3], s_n[3], s_t[3];
  double u[3], u1[3], u2[3], u3[3];
  double l, l2, r, rmin, rmax, f2;
  double E, F, bt;
  double Tl, Imin, Imax, f3;
  double d1, d2;

  // These are optimized curvatures suitable for the human brain
  rmin = self->GetRMin();
  rmax = self->GetRMax();

  // These are used to compute the smoothness update fraction
  E = 0.5*(1/rmin + 1/rmax);
  F = 6.0/(1/rmin - 1/rmax);

  // BT
  bt = self->GetBT();

  // Min/Max distances for the search
  d1 = self->GetD1();
  d2 = self->GetD2();

  // Find and store neighbours

  // Switching from VTK containers to STL containers
  // because they're much faster for iterating over

  // Points
  myPoint pt;
  std::vector<myPoint> brainPoints(nPoints);

  // Neighbouring Points Ids
  std::vector<int> thisPointNeighbourIds; // 5 or 6
  std::vector< std::vector<int> > pointNeighbourIds(nPoints);
  std::vector<int>::iterator neighbourIdIter;

  // Neighbouring Cells (their vertex point ids)
  myVertIds vIds;
  std::vector<myVertIds> thisPointNeighbourVertIds; // 5 or 6
  std::vector< std::vector<myVertIds> > pointNeighbourVertIds(nPoints);
  std::vector<myVertIds>::iterator neighbourVertIdsIter;

  // Loop over each point in brainPolyData and save it into an STL vector.
  // Also, use vtkPolyData functions to figure out where our neighbours are
  for (vtkIdType ptId = 0; ptId < nPoints ; ptId++)
    {
    // Our target point
    double point[3];
    brainPolyData->GetPoint(ptId, point);

    pt = point; // pass the c array into a c++ class for the vector
    brainPoints[ptId] = pt;

    // The cells attached to the target point
    unsigned short nCells;
    vtkIdType *pointCells;
    brainPolyData->GetPointCells( ptId, nCells, pointCells );

    // The points inside the attached cells
    for (unsigned short cellIdx = 0; cellIdx < nCells; cellIdx++)
      {
      vtkIdType cellId = pointCells[cellIdx];
      vtkIdType npts, *pts;
      brainPolyData->GetCellPoints(cellId, npts, pts);

      // If the point is not our target point, it's a neighbour
      for (vtkIdType nIdx = 0; nIdx < npts ; nIdx++ )
        {
        vtkIdType nId = pts[nIdx];
        vIds[nIdx] = nId;

        if (ptId != nId)
          {
          if (std::find( thisPointNeighbourIds.begin(),
                         thisPointNeighbourIds.end(),
                         nId) == thisPointNeighbourIds.end()) // not found
            {
            thisPointNeighbourIds.push_back(nId);
            }
          }
        }
      thisPointNeighbourVertIds.push_back(vIds);
      }
    pointNeighbourIds[ptId] = thisPointNeighbourIds;
    thisPointNeighbourIds.clear();

    pointNeighbourVertIds[ptId] = thisPointNeighbourVertIds;
    thisPointNeighbourVertIds.clear();
    }

  // Temp variables
  myPoint target, neighbour;
  myPoint v0, v1, v2;
  double temp0[3], temp1[3], temp2[3];
  double cross[3];
  double incX, incY, incZ;
  double u3multiplier;

  // avoid compiler warnings
  Tl = l = l2 = Imin = Imax = 0.0;
  incX = incY = incZ = 0.0;
  nNeighbours = 0;
  u3multiplier = 1.0;

  // We need somewhere to put the new points during each iteration
  std::vector<myPoint> updatePoints(nPoints);

  // THE MAIN LOOP
  vtkIdType idx;
  vtkIdType incs[3];
  inData->GetIncrements(incs);
  iteration = 0;
  nIterations = self->GetNumberOfIterations();

  std::vector<myPoint>::iterator ptIter;
  std::vector< std::vector<int> >::iterator nIter;
  std::vector< std::vector<myVertIds> >::iterator vIter;
  std::vector<myPoint>::iterator uIter;

  // Copy brain points as original points
  vtkPoints *originalPoints = vtkPoints::New();
  for (ptIter = brainPoints.begin();
       ptIter != brainPoints.end();
       ptIter++)
    {
    target = *ptIter;
    originalPoints->InsertNextPoint(target.xyz);
    }

  while (iteration < nIterations)
    {
    // Update l every 50 iterations
    if (iteration%50 == 0)
      {
      suml2 = 0.0;

      nIter = pointNeighbourIds.begin();
      for (ptIter = brainPoints.begin();
           ptIter != brainPoints.end();
           ptIter++, nIter++)
        {
        target = *ptIter;
        thisPointNeighbourIds = *nIter;

        nNeighbours = static_cast<int>(thisPointNeighbourIds.size());

        this_suml2 = 0.0;
        for (neighbourIdIter = thisPointNeighbourIds.begin();
             neighbourIdIter!=thisPointNeighbourIds.end(); neighbourIdIter++)
          {
          neighbour = brainPoints[*neighbourIdIter];
          this_suml2 += vtkMath::Distance2BetweenPoints(target.xyz,
                                                        neighbour.xyz);
          }

        suml2 += this_suml2/nNeighbours;
        }
      l2 = suml2/static_cast<double>(nPoints);
      l = sqrt(l2);
      u3multiplier = 0.05*l; // just calc this when l changes
      }

    nIter = pointNeighbourIds.begin();
    vIter = pointNeighbourVertIds.begin();
    uIter = updatePoints.begin();
    for (ptIter = brainPoints.begin();
         ptIter != brainPoints.end();
         ptIter++, nIter++, vIter++, uIter++)
      {
      target = *ptIter;
      thisPointNeighbourIds = *nIter;
      thisPointNeighbourVertIds = *vIter;

      // reset the normal and average location of our neighbours
      n_hat[0] = n_hat[1] = n_hat[2] = 0.0;
      temp2[0] = temp2[1] = temp2[2] = 0.0;
      nNeighbours = static_cast<int>(thisPointNeighbourIds.size());

      //Normal
      // loop over each neighour's vertex ids
      for( neighbourVertIdsIter =  thisPointNeighbourVertIds.begin();
           neighbourVertIdsIter != thisPointNeighbourVertIds.end();
           neighbourVertIdsIter++)
        {
        // 3 Ids that give us the location of our neighbours verts
        vIds = *neighbourVertIdsIter;

        v0 = brainPoints[vIds[0]];
        v1 = brainPoints[vIds[1]];
        v2 = brainPoints[vIds[2]];

        // Find the normal of our triangle.  We used vtkPolyData->GetPointCells
        // to find our neighbours so we can be vtk-certain that v1
        // is our target
        temp0[0] = v2.xyz[0] - v1.xyz[0];
        temp0[1] = v2.xyz[1] - v1.xyz[1];
        temp0[2] = v2.xyz[2] - v1.xyz[2];

        temp1[0] = v0.xyz[0] - v1.xyz[0];
        temp1[1] = v0.xyz[1] - v1.xyz[1];
        temp1[2] = v0.xyz[2] - v1.xyz[2];

        vtkMath::Cross(temp0, temp1, cross);

        n_hat[0] += cross[0];
        n_hat[1] += cross[1];
        n_hat[2] += cross[2];
        }

      vtkMath::Normalize(n_hat);

      // Mean location
      for (neighbourIdIter = thisPointNeighbourIds.begin();
           neighbourIdIter != thisPointNeighbourIds.end();
           neighbourIdIter++)
        {
        neighbour = brainPoints[*neighbourIdIter];
        temp2[0] += neighbour.xyz[0];
        temp2[1] += neighbour.xyz[1];
        temp2[2] += neighbour.xyz[2];
        }

      // s is a vector that takes target to the mean position of
      // it's neighbours
      s[0] = temp2[0]/nNeighbours - target.xyz[0];
      s[1] = temp2[1]/nNeighbours - target.xyz[1];
      s[2] = temp2[2]/nNeighbours - target.xyz[2];

      // s is decomposed into a normal component s_n
      // and a tangential component s_t
      dotp = vtkMath::Dot(s,n_hat);
      s_n[0] = dotp*n_hat[0];
      s_n[1] = dotp*n_hat[1];
      s_n[2] = dotp*n_hat[2];

      s_t[0] = s[0]-s_n[0];
      s_t[1] = s[1]-s_n[1];
      s_t[2] = s[2]-s_n[2];

      // Target is moved by update vector u which is made up of u1+u2+u3

      // Update component 1: within-surface vertex spacing
      u1[0] = 0.5*s_t[0];
      u1[1] = 0.5*s_t[1];
      u1[2] = 0.5*s_t[2];

      // Update component 2: surface smothness control

      // the sigmoid function
      //r = l2/(2*vtkMath::Norm(s_n));
      r = 2*vtkMath::Norm(s_n)/l2; // actually the inverse

      // the smoothness fraction
      //f2 = 0.5*(1+tanh(F*(1/r-E)));
      f2 = 0.5*(1+tanh(F*(r-E))); // use r inverse

      // u2
      u2[0] = f2*s_n[0];
      u2[1] = f2*s_n[1];
      u2[2] = f2*s_n[2];

      // Update component 3: interact with the image data
      u3[0] = u3[1] = u3[2] = 0.0;
      f3 = 0.0;

      // Direction is anti-parallel to normal
      direction[0] = -n_hat[0];
      direction[1] = -n_hat[1];
      direction[2] = -n_hat[2];

      // Start the search 1mm in from target
      fstart[0] = target.xyz[0] + direction[0];
      fstart[1] = target.xyz[1] + direction[1];
      fstart[2] = target.xyz[2] + direction[2];

      // Get the far end of the search
      fend[0] = fstart[0] + (d1-1.0)*direction[0];
      fend[1] = fstart[1] + (d1-1.0)*direction[1];
      fend[2] = fstart[2] + (d1-1.0)*direction[2];

      //Change to voxel coordinates
      location[0] = (fstart[0]-origin[0])/spacing[0];
      location[1] = (fstart[1]-origin[1])/spacing[1];
      location[2] = (fstart[2]-origin[2])/spacing[2];

      start[0] = static_cast<int>(location[0] + 0.5);
      start[1] = static_cast<int>(location[1] + 0.5);
      start[2] = static_cast<int>(location[2] + 0.5);

      end[0] = static_cast<int>((fend[0]-origin[0])/spacing[0] + 0.5);
      end[1] = static_cast<int>((fend[1]-origin[1])/spacing[1] + 0.5);
      end[2] = static_cast<int>((fend[2]-origin[2])/spacing[2] + 0.5);

      // If the search remains inside the volume, continue.
      if (extent[0] <= start[0] && start[0] <= extent[1] &&
          extent[2] <= start[1] && start[1] <= extent[3] &&
          extent[4] <= start[2] && start[2] <= extent[5] &&
          extent[0] <= end[0]   && end[0] <= extent[1] &&
          extent[2] <= end[1]   && end[1] <= extent[3] &&
          extent[4] <= end[2]   && end[2] <= extent[5])
        {
        Imin = Tm;
        Imax = TH;

        idx = (start[0]*incs[0]+
               start[1]*incs[1]+
               start[2]*incs[2]);

        value = static_cast<double>(inPtr[idx]);

        Imin = std::min(Imin, value);
        Imax = std::max(Imax, value);

        incX = direction[0]/spacing[0];
        incY = direction[1]/spacing[1];
        incZ = direction[2]/spacing[2];

        d = 1.0;
        // this loop causes alot of stalling on G4 PPCs.
        while (d < d1)
          {
          location[0] += incX;
          location[1] += incY;
          location[2] += incZ;

          idx = (static_cast<int>(location[0] + 0.5)*incs[0]+
                 static_cast<int>(location[1] + 0.5)*incs[1]+
                 static_cast<int>(location[2] + 0.5)*incs[2]);
          value = static_cast<double>(inPtr[idx]); // here's the stall

          // Min search up to d1
          Imin = std::min(Imin, value);

          // Max search up to d2
          if (d<d2)
            {
            Imax = std::max(Imax, value);
            }

          d += 1.0; // Step by 1mm?
          }

        Imin = std::max(T2, Imin);
        Imax = std::min(Tm, Imax);

        // The local background
        Tl = (Imax-T2)*bt+T2;

        // Decide which way u3 should go
        if ((Imax-T2) > 0)
          {
          f3 = 2.0*(Imin-Tl)/(Imax-T2);
          }
        else
          {
          f3 = 2.0*(Imin-Tl);
          }

        // NOTE: in the paper (eq. 13) this is s_n_hat, the wrong direction!
        u3[0] = u3multiplier*f3*n_hat[0];
        u3[1] = u3multiplier*f3*n_hat[1];
        u3[2] = u3multiplier*f3*n_hat[2];
        }

      //  The final update vector
      u[0] = u1[0] + u2[0] + u3[0];
      u[1] = u1[1] + u2[1] + u3[1];
      u[2] = u1[2] + u2[2] + u3[2];

      // New target
      temp0[0] = target.xyz[0] + u[0];
      temp0[1] = target.xyz[1] + u[1];
      temp0[2] = target.xyz[2] + u[2];

      *uIter = temp0;
      }
    // we've got new points
    std::copy(updatePoints.begin(), updatePoints.end(), brainPoints.begin());
    iteration++;
    }
  // Switch back to VTK containers
  vtkPoints *newPoints = brainPolyData->GetPoints();
  for (vtkIdType ptId = 0; ptId < nPoints; ptId++)
    {
    target = brainPoints[ptId];
    newPoints->SetPoint( ptId, target.xyz );
    }

  brainPolyData->Modified();

  // Aviod ugly poly data - unnecessary?
  vtkCleanPolyData *cleanPoly = vtkCleanPolyData::New();
  cleanPoly->SET_INPUT_DATA(brainPolyData);
  cleanPoly->Update();

  self->GetBrainMesh()->ShallowCopy(cleanPoly->GetOutput());

  //Use the brain mesh to stencil out the non-brain
  vtkPolyDataToImageStencil *theStencil = vtkPolyDataToImageStencil::New();
  vtkImageStencil *imageStencil = vtkImageStencil::New();

  theStencil->SET_INPUT_DATA(self->GetBrainMesh());
  theStencil->SetInformationInput(inData);
  theStencil->Update();

  imageStencil->SET_STENCIL_DATA(theStencil->GetOutput());
  imageStencil->SET_INPUT_DATA(inData);
  imageStencil->SetBackgroundValue(T2);
  imageStencil->Update();

  self->GetOutput()->ShallowCopy(imageStencil->GetOutput());

  //Clean up
  originalPoints->Delete();
  theStencil->Delete();
  imageStencil->Delete();
  brainPolyData->Delete();
  pointNeighbourList->Delete();
  cleanPoly->Delete();
}

} // end anonymous namespace

//----------------------------------------------------------------------------
int vtkImageMRIBrainExtractor::RequestData(
  vtkInformation *, vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkImageData *inData = vtkImageData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkImageData *outData = vtkImageData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  int outExt[6];
  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), outExt);

#if VTK_MAJOR_VERSION >= 6
  this->AllocateOutputData(outData, outInfo, outExt);
#else
  this->AllocateOutputData(outData, outExt);
#endif

  
  void *inPtr = inData->GetScalarPointerForExtent(outExt);

  if (inData->GetScalarType() != outData->GetScalarType())
    {
    vtkErrorMacro("Execute: Output ScalarType "
                  << outData->GetScalarType()
                  << ", must Input ScalarType "
                  << inData->GetScalarType());
    return 0;
    }

  switch (inData->GetScalarType())
    {
    case VTK_UNSIGNED_CHAR:
      vtkImageMRIBrainExtractorExecute(
        this, inData, static_cast<unsigned char *>(inPtr));
      break;
    case VTK_SHORT:
      vtkImageMRIBrainExtractorExecute(
        this, inData, static_cast<short *>(inPtr));
      break;
    case VTK_UNSIGNED_SHORT:
      vtkImageMRIBrainExtractorExecute(
        this, inData, static_cast<unsigned short *>(inPtr));
      break;
    default:
      vtkErrorMacro(<< "Execute: "
        "Requires short, unsigned short, or unsigned char");
      return 0;
    }

  return 1;
}

//----------------------------------------------------------------------------
void vtkImageMRIBrainExtractor::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "BrainExtent: "
     << this->BrainExtent[0] << " " << this->BrainExtent[1] << " "
     << this->BrainExtent[2] << " " << this->BrainExtent[3] << " "
     << this->BrainExtent[4] << " " << this->BrainExtent[5] << "\n";
  os << indent << "BT: " << this->BT << "\n";
  os << indent << "D1: " << this->D1 << "\n";
  os << indent << "D2: " << this->D2 << "\n";
  os << indent << "RMin: " << this->RMin << "\n";
  os << indent << "RMax: " << this->RMax << "\n";
  os << indent << "NumberOfIterations: "
     << this->NumberOfIterations << "\n";
  os << indent << "NumberOfTessellations: "
     << this->NumberOfTessellations << "\n";
}
