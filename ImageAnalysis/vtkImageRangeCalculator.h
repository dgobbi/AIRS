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
// .Name vtkImageRangeCalculator - Calculate Data Range.
// .SECTION Description
// vtkImageRangeCalculator calculates the Data Range for a given
// fraction range.  i.e. given a fraction range of (0.05, 0.95), this
// class will return the data values that correspond to 5% and 95% of
// the total histogram area.

#ifndef __vtkImageRangeCalculator_h
#define __vtkImageRangeCalculator_h

#include "vtkSystemIncludes.h"
#ifndef vtkFloatingPointType
#define vtkFloatingPointType vtkFloatingPointType
typedef float vtkFloatingPointType;
#endif

#include "vtkImageToImageFilter.h"

class vtkImageStencilData;

class VTK_EXPORT vtkImageRangeCalculator : public vtkImageToImageFilter
{
public:
  static vtkImageRangeCalculator* New();
  vtkTypeMacro(vtkImageRangeCalculator, vtkImageToImageFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  vtkSetVector2Macro(AreaFractionRange, double);
  vtkGetVector2Macro(AreaFractionRange, double);

  vtkGetVector2Macro(DataRange, double);

  void SetNumberOfThreads(int i);

  void Calculate();

protected:
  vtkImageRangeCalculator();
  ~vtkImageRangeCalculator();

  void ThreadedExecute(vtkImageData *inData, vtkImageData *outData, 
                       int ext[6], int id);

  double AreaFractionRange[2];
  double DataRange[2];

private:
  vtkImageRangeCalculator(const vtkImageRangeCalculator&);  // Not implemented.
  void operator=(const vtkImageRangeCalculator&);  // Not implemented.
};

#endif
