import sys, os

import vtk

if os.name == 'posix':
    import libvtkAtamaiImagingPython
    vtk.vtkImageProjection = libvtkAtamaiImagingPython.vtkImageProjection
else:
    import vtkAtamaiImagingPython
    vtk.vtkImageProjection = vtkAtamaiImagingPython.vtkImageProjection

project = vtk.vtkImageProjection()

sys.exit(0)

