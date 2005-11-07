import sys, os

import vtk
import libvtkAtamaiImagingPython
vtk.vtkImageProjection = libvtkAtamaiImagingPython.vtkImageProjection

project = vtk.vtkImageProjection()

sys.exit(0)

