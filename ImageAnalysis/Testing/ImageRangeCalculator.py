#!/usr/bin/env python

import os, random, string
import vtk
from vtk.test import Testing

if os.name == 'posix':
    import libvtkAtamaiImagingPython
    vtk.vtkImageRangeCalculator = libvtkAtamaiImagingPython.vtkImageRangeCalculator
else:
    import vtkAtamaiImagingPython
    vtk.vtkImageRangeCalculator = vtkAtamaiImagingPython.vtkImageRangeCalculator

class ImageRangeCalculate(Testing.vtkTest):
    def __init__(self, test):
        Testing.vtkTest.__init__(self, test)

    def test(self):

        # setup some test data
        source = vtk.vtkImageGaussianSource()
        source.SetWholeExtent(0, 63, 0, 63, 0, 63)
        source.SetCenter(31.5, 31.5, 31.5)
        source.SetMaximum(255)
        source.SetStandardDeviation(32)

        # cast to unsigned char
        cast = vtk.vtkImageShiftScale()
        cast.SetInput(source.GetOutput())
        cast.ClampOverflowOn()
        cast.SetShift(-128.0)
        cast.SetScale(1.0)
        cast.SetOutputScalarTypeToDouble()

        imageData = cast.GetOutput()

        # add a stencil
        sphere = vtk.vtkSphere()
        sphere.SetCenter(31.5, 31.5, 31.5)
        sphere.SetRadius(32)

        functionToStencil = vtk.vtkImplicitFunctionToImageStencil()
        functionToStencil.SetInput(sphere)
        
        functionToStencil.Update() # REQUIRED!

        stencilData = functionToStencil.GetOutput()

        # accumulate the histogram
        accumulate = vtk.vtkImageAccumulate()
        accumulate.SetInput(imageData)
        accumulate.SetStencil(stencilData)
        accumulate.SetComponentSpacing(1, 0, 0)
        accumulate.SetComponentOrigin(-128.0, 0, 0)
        accumulate.SetComponentExtent(0, 255, 0, 0, 0, 0)

        histogramData = accumulate.GetOutput()

        # calculate values for a range
        range = vtk.vtkImageRangeCalculator()
        range.SetInput(histogramData)
        range.SetAreaFractionRange(0.05, 0.95)

        # This is the same:
        # rangeData = range.GetOutput()
        # rangeData.Update()
        # as this:
        range.Calculate()

        # Get the data range and check it
        dataRange = range.GetDataRange()
        
        correctResult = (dataRange[0] == 28.0 and dataRange[1] == 109.0)
        self.failUnless( correctResult, "Data range is incorrect")
        
        
if __name__ == "__main__":
    Testing.main([(ImageRangeCalculate, 'test')])
