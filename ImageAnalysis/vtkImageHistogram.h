/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageHistogram.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageHistogram - Generate the histogram for an image.
// .SECTION Description
// vtkImageHistogram generates a histogram from its input, and provides
// a 2D image of the histogram as its output.  Unlike vtkImageAccumulate,
// a multi-component image does not result in a multi-dimensional
// histogram.  Instead, the histogram is computed over all components.

#ifndef __vtkImageHistogram_h
#define __vtkImageHistogram_h

#include "vtkThreadedImageAlgorithm.h"

class vtkImageStencilData;
class vtkIdTypeArray;

class VTK_EXPORT vtkImageHistogram : public vtkThreadedImageAlgorithm
{
public:
  static vtkImageHistogram *New();
  vtkTypeMacro(vtkImageHistogram,vtkThreadedImageAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // If this is On, then the histogram binning will be done automatically.
  // For char and unsigned char data, there will be 256 bins with unit
  // spacing.  For data of type short and larger, there will be between
  // 256 and MaximumNumberOfBins, depending on the range of the data, and
  // the BinOrigin will be set to zero if no negative values are present,
  // or to the smallest negative value if negative values are present.
  // For float data, the MaximumNumberOfBins will always be used.
  // The BinOrigin and BinSpacing will be set so that they provide a mapping
  // from bin index to scalar value.
  vtkSetMacro(AutomaticBinning, int);
  vtkBooleanMacro(AutomaticBinning, int);
  vtkGetMacro(AutomaticBinning, int);

  // Description:
  // The maximum number of bins to use when AutomaticBinning is On.
  // The default value is 65536, enough to capture the full range of
  // short and unsigned short values.
  vtkSetMacro(MaximumNumberOfBins, int);
  vtkGetMacro(MaximumNumberOfBins, int);

  // Description:
  // The number of bins in histogram (default 256).  This is automatically
  // computed if AutomaticBinning is On. 
  vtkSetMacro(NumberOfBins, int);
  vtkGetMacro(NumberOfBins, int);

  // Description:
  // The value for the center of the first bin (default 0).  This is
  // automatically computed if AutomaticBinning is On.
  vtkSetMacro(BinOrigin, double);
  vtkGetMacro(BinOrigin, double);

  // Description:
  // The bin spacing (default 1).  This is automatically computed if
  // AutomaticBinning is On.
  vtkSetMacro(BinSpacing, double);
  vtkGetMacro(BinSpacing, double);

  // Description:
  // Use a stencil to compute the histogram for just a part of the image.
  void SetStencil(vtkImageStencilData *stencil);
  vtkImageStencilData *GetStencil();

  // Description:
  // If this is On, then a histogram image will be produced as the output.
  // Otherwise, the histogram is only returned as a vtkDataArray.
  vtkSetMacro(GenerateHistogramImage, int);
  vtkBooleanMacro(GenerateHistogramImage, int);
  vtkGetMacro(GenerateHistogramImage, int);

  // Description:
  // Set the size of the histogram image that is produced as output.
  // The default is 256 by 256.
  vtkSetVector2Macro(HistogramImageSize, int);
  vtkGetVector2Macro(HistogramImageSize, int);

  // Description:
  // Get the histogram as a vtkIdTypeArray.  You must call Update()
  // before calling this method.
  vtkDataArray *GetHistogram();

  // Description:
  // Get the total count of the histogram.  This will be the number of
  // voxels times the number of components.
  vtkIdType GetTotal() { return this->Total; }

  // Description:
  // This is part of the executive, but is public so that it can be accessed
  // by non-member functions.
  virtual void ThreadedRequestData(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector,
                                   vtkImageData ***inData,
                                   vtkImageData **outData, int ext[6], int id);
protected:
  vtkImageHistogram();
  ~vtkImageHistogram();

  virtual int RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inInfo,
                                 vtkInformationVector *vtkNotUsed(outInfo));
  virtual int RequestInformation(vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inInfo,
                                 vtkInformationVector *vtkNotUsed(outInfo));
  virtual int RequestData(vtkInformation *,
			  vtkInformationVector **,
			  vtkInformationVector *);

  virtual int FillInputPortInformation(int port, vtkInformation *info);
  virtual int FillOutputPortInformation(int port, vtkInformation *info);

  // Description:
  // Compute the range of the data.  The GetScalarRange() function of
  // vtkImageData only computes the range of the first component, but
  // this filter requires the range for all components.
  void ComputeImageScalarRange(vtkImageData *data, double range[2]);

  // Description:
  // This is a hook function for subclasses that analyze the histogram.
  virtual void ComputeStatistics();

  int AutomaticBinning;
  int MaximumNumberOfBins;
  
  int HistogramImageSize[2];
  int GenerateHistogramImage;

  int NumberOfBins;
  double BinOrigin;
  double BinSpacing;

  vtkIdTypeArray *Histogram;
  vtkIdType Total;

  vtkIdType *ThreadOutput[VTK_MAX_THREADS];
  int ThreadBinRange[VTK_MAX_THREADS][2];

private:
  vtkImageHistogram(const vtkImageHistogram&);  // Not implemented.
  void operator=(const vtkImageHistogram&);  // Not implemented.
};

#endif
