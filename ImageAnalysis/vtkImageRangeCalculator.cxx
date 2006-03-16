/*=========================================================================

Copyright (c) 2006 Atamai, Inc.

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
#include "vtkImageRangeCalculator.h"

#include "vtkObjectFactory.h"
#include "vtkImageStencilData.h"
#include "vtkImageData.h"

//vtkCxxRevisionMacro(vtkImageRangeCalculator, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkImageRangeCalculator);

//----------------------------------------------------------------------------
vtkImageRangeCalculator::vtkImageRangeCalculator()
{
  // no threads
  this->Superclass::SetNumberOfThreads(1);

  // ivars
  this->AreaFractionRange[0] = 0.0;
  this->AreaFractionRange[1] = 1.0;
  this->DataRange[0] = 0; // ??
  this->DataRange[1] = 0; // ??
  
}

//----------------------------------------------------------------------------
vtkImageRangeCalculator::~vtkImageRangeCalculator()
{
}

//----------------------------------------------------------------------------
void vtkImageRangeCalculator::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
// Override and ignore
void vtkImageRangeCalculator::SetNumberOfThreads(int i)
{
  return;
}


//----------------------------------------------------------------------------
double* vtkImageRangeCalculatorExecute(vtkImageRangeCalculator* self, 
				      vtkImageData *inData, void *inPtr,
				      int inExt[6], double Range[2])
{
  // TODO:
  // - Add support for 2D and 3D histograms
  // - Fix the progress methods

  int numscalars;
  int idX, idY, idZ;
  int idXmin, idXmax, iter;
  vtkIdType inIncX, inIncY, inIncZ;
  int scalarSize;
  vtkIdType inInc[3];
  unsigned long count = 0;
  unsigned long target;
  double point[4];
  double f;
  vtkFloatingPointType *inSpacing, *inOrigin;

  // for conversion to data coordinates
  inOrigin = inData->GetOrigin();
  inSpacing = inData->GetSpacing();
  
  // for the progress meter
  target = (unsigned long)
    ((inExt[5]-inExt[4]+1)*(inExt[3]-inExt[2]+1)/50.0);
  target++;
  
  // Get Increments to march through data 
  int tmpInc[3]; // in VTK 4.4 and earlier, increments are int
  inData->GetIncrements(tmpInc);
  inInc[0] = tmpInc[0];
  inInc[1] = tmpInc[1];
  inInc[2] = tmpInc[2];
  scalarSize = inData->GetScalarSize();
  numscalars = inData->GetNumberOfScalarComponents();

  // The output of vtkImageAccumulate is always int
  int *tmpPtr = (int *)inPtr;

  // Loop 1: Find the total area
  int *tmpPtr1 = tmpPtr;
  int nbins = inExt[1]-inExt[2]; // ASSUMING 1D
  int area = 0;
  for (idZ = inExt[4]; idZ <= inExt[5]; idZ++)
    {
    for (idY = inExt[2]; idY <= inExt[3]; idY++)
      {
      // update the progress if this is the main thread
      if (!(count%target)) 
	{
	self->UpdateProgress(count/(50.0*target));
	}
      count++;

      for (idX = inExt[0]; idX <= inExt[1]; idX++)
	{
	area += *tmpPtr1;
	tmpPtr1++;
	} 
      }
    }

  double *fractionRange = self->GetAreaFractionRange();

  int minArea, maxArea;
  minArea = int(fractionRange[0]*area+0.5);
  maxArea = int(fractionRange[1]*area+0.5);

  int tempArea, minFound, maxFound;
  tempArea = minFound = maxFound = 0;

  int minBin, maxBin;
  minBin = maxBin = 0;

  // Loop 2: Find the data range
  int *tmpPtr2 = tmpPtr;
  for (idZ = inExt[4]; idZ <= inExt[5]; idZ++)
    {
    for (idY = inExt[2]; idY <= inExt[3]; idY++)
      {
      // update the progress if this is the main thread
      if (!(count%target)) 
	{
	self->UpdateProgress(count/(50.0*target));
	}
      count++;

      for (idX = inExt[0]; idX <= inExt[1]; idX++)
	{
	tempArea += *tmpPtr2;

	// new min
	if (tempArea <= minArea)
	  {
	  minBin = idX;
	  }

	// new max
	if (tempArea <= maxArea)
	  {
	  maxBin = idX;
	  }

	tmpPtr2++;
	} 
      }
    }
 
  Range[0] = (minBin+inOrigin[0])*inSpacing[0];
  Range[1] = (maxBin+inOrigin[0])*inSpacing[0];

}

// ---------------------------------------------------------------------------
// This method is passed a input and output region, and executes the filter
// algorithm to fill the output from the input.
// It just executes a switch statement to call the correct function for
// the regions data types.
void vtkImageRangeCalculator::ThreadedExecute(vtkImageData *inData, 
                                      vtkImageData *outData,
                                      int outExt[6], int id)
{
  int inExt[6];
  inData->GetExtent(inExt);
  void *inPtr = outData->GetScalarPointerForExtent(inExt);
  if (inExt[0] <= inExt[1] && inExt[2] <= inExt[3] && inExt[4] <= inExt[5])
    {
    inPtr = inData->GetScalarPointerForExtent(inExt);
    }
  
  vtkDebugMacro(<< "Execute: inData = " << inData 
                << ", outData = " << outData);

  vtkImageRangeCalculatorExecute(this, inData, inPtr, inExt,
				 this->DataRange);
  
}

void vtkImageRangeCalculator::Calculate()
{
  // Call ThreadedExecute
  vtkImageData *outData = this->GetOutput();
  outData->UpdateInformation();
  outData->Update();
}

