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
// .NAME vtkImageMRIBrainExtractor - Create brain mask for an MR image.
// .SECTION Description
// vtkImageMRIBrainExtractor uses a heuristic approach(1) to segment brain
// from non-brain in T1-weighted MRI images.  Use the GetOutput method to get
// the segmented (brain only) vtkImageData and the GetBrainMesh to get the
// vtkPolyData outlining the brain surface.
//
// 1. Smith, S.M., "Fast Robust Automated Brain Extraction,"
// Human Brain Mapping, 17:143-155, 2002.


#ifndef vtkImageMRIBrainExtractor_h
#define vtkImageMRIBrainExtractor_h

#include "vtkImageSegmentationModule.h" // For export macro
#include "vtkImageAlgorithm.h"

class vtkPolyData;
class vtkPoints;

class VTKIMAGESEGMENTATION_EXPORT vtkImageMRIBrainExtractor :
  public vtkImageAlgorithm
{
public:
  static vtkImageMRIBrainExtractor* New();
  vtkTypeMacro(vtkImageMRIBrainExtractor, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the region where the brain is expected to be.  If this is
  // not set, or if it is set to be larger than the full extent of
  // the image, then the full extent of the image will be used.  
  vtkSetVector6Macro(BrainExtent, int);
  vtkGetVector6Macro(BrainExtent, int);

  // Description:
  // Get the vtkPolyData representing the brain surface
  vtkPolyData *GetBrainMesh();

  // Description:
  // BT (0.0, 1.0) is a fractional constant that controls the fit of the
  // brain surface segmentation.  0.5 is the default.
  vtkSetMacro(BT, double);
  vtkGetMacro(BT, double);

  // Description:
  // D1 is the distance (mm) that is searched for minimum intensity values
  // during segmentation.  Typically D1 = 2*D2. Default is 20mm.
  vtkSetMacro(D1, double);
  vtkGetMacro(D1, double);

  // Description:
  // D2 is the distance (mm) that is searched for maximum intensity values
  // during segmentation.  Typically D2 = 0.5*D1. Default is 10mm.
  vtkSetMacro(D2, double);
  vtkGetMacro(D2, double);

  // Description:
  // Number of iterations used to deform the brain mesh. Default is 1000.
  vtkSetMacro(NumberOfIterations, int);
  vtkGetMacro(NumberOfIterations, int);

  // Description:
  // Number of tessellations of an icosahedron (20-sided polygon).  The
  // default is 4, which gives 2562 points over the surface of the brain.
  vtkSetMacro(NumberOfTessellations, int);
  vtkGetMacro(NumberOfTessellations, int);

  // Description:
  // The minimum acceptable radius of curvature.  A local radius of
  // curvature below this value will cause increased smoothing. The
  // default is 3.0.
  vtkSetMacro(RMin, double);
  vtkGetMacro(RMin, double);

  // Description:
  // The maximum acceptable radius of curvature.  A local radius of
  // curvature above RMin and below RMax is smoothed according to a
  // non-linear function that decreases from RMin to RMax.  A local
  // radius of curvature about RMax will result in no smoothing.
  // The default is 10.0.
  vtkSetMacro(RMax, double);
  vtkGetMacro(RMax, double);

protected:
  vtkImageMRIBrainExtractor();
  ~vtkImageMRIBrainExtractor();

  virtual int RequestData(vtkInformation *,
			  vtkInformationVector **,
			  vtkInformationVector *);

  static void ComputeMeshCentroid(vtkPolyData *data, double cen[3]);

  int BrainExtent[6];
  vtkPolyData *BrainMesh;
  double BT;
  int NumberOfIterations;
  int NumberOfTessellations;
  double D1;
  double D2;
  double RMin;
  double RMax;

private:
  vtkImageMRIBrainExtractor(const vtkImageMRIBrainExtractor&);  // Not implemented.
  void operator=(const vtkImageMRIBrainExtractor&);  // Not implemented.
};

#endif
