#!/usr/bin/env python

import os, random, string
import vtk
from vtk.test import Testing

if os.name == 'posix':
    import libvtkAtamaiRegistrationPython
    vtk.vtkImageRegistration = libvtkAtamaiRegistrationPython.vtkImageRegistration
else:
    import vtkAtamaiRegistrationPython
    vtk.vtkImageRegistration = vtkAtamaiRegistrationPython.vtkImageRegistration

# constants for the optimizers
VTK_OPTIMIZER_AMOEBA                           = 0

# constants for the metrics
VTK_METRIC_NORMALIZED_CROSS_CORRELATION        = 0
VTK_METRIC_NORMALIZED_MUTUAL_INFORMATION       = 1

# constants for the interpolators
VTK_INTERPOLATOR_NEAREST_NEIGHBOR              = 0
VTK_INTERPOLATOR_LINEAR                        = 1

# constants for the transform types
VTK_TRANSFORM_CENTERED                         = 0

TRANSLATION_X                                  = 3.0
TRANSLATION_Y                                  = 5.0
TRANSLATION_Z                                  = 7.0

class TestImageRegistration(Testing.vtkTest):
    def __init__(self, test):
        """
        """
        Testing.vtkTest.__init__(self, test)

    def test(self):
        """
        """
        #setup test data
        source = vtk.vtkImageGaussianSource()
        source.SetWholeExtent(0, 64, 0, 64, 0, 64)
        source.SetCenter(31.5, 31.5, 31.5)
        source.SetMaximum(255)
        source.SetStandardDeviation(32)

        sourceCaster = vtk.vtkImageShiftScale()
        sourceCaster.SetInput( source.GetOutput() )
        sourceCaster.SetOutputScalarTypeToFloat()
        sourceCaster.UpdateWholeExtent()
        
        sourceImage = sourceCaster.GetOutput()

        transform = vtk.vtkTransform()
        transform.Translate(TRANSLATION_X, TRANSLATION_Y, TRANSLATION_Z)
        
        reslicer = vtk.vtkImageReslice()
        reslicer.SetInput(sourceImage)
        reslicer.SetResliceTransform(transform)

        targetImage = reslicer.GetOutput()

        registration = vtk.vtkImageRegistration()
        registration.SetFixedImage(targetImage)
        registration.SetMovingImage(sourceImage)

        registration.SetOptimizerType( VTK_OPTIMIZER_AMOEBA )
        registration.SetOptimizerParameter( "MaxNumberOfIterations", 100 )

        registration.SetMetricType( VTK_METRIC_NORMALIZED_MUTUAL_INFORMATION )
        registration.SetMetricParameter( "AComponentExtent", 32 )
        registration.SetMetricParameter( "BComponentExtent", 32 )
        
        registration.SetTransformType( VTK_TRANSFORM_CENTERED )
        registration.SetTransformParameter( "CenterX", 31.5 )
        registration.SetTransformParameter( "CenterY", 31.5 )
        registration.SetTransformParameter( "CenterZ", 31.5 )
        registration.SetTransformParameter( "TranslationX", 0.0 )
        registration.SetTransformParameter( "TranslationY", 0.0 )
        registration.SetTransformParameter( "TranslationZ", 0.0 )
        registration.SetTransformParameter( "RotationAngleX", 0.0 )
        registration.SetTransformParameter( "RotationAngleY", 0.0 )
        registration.SetTransformParameter( "RotationAngleZ", 0.0 )
        registration.SetTransformParameter( "IsotropicScale", 1.0 )
        registration.SetTransformParameter( "TranslationScale", 10.0 )
        registration.SetTransformParameter( "RotationScale", 0.0 )
        registration.SetTransformParameter( "ScalingScale", 0.0 )
        
        registration.SetInterpolatorType( VTK_INTERPOLATOR_LINEAR )
        
        registration.UpdateRegistration()

        tx = registration.GetTransformParameter( "TranslationX" )
        ty = registration.GetTransformParameter( "TranslationY" )
        tz = registration.GetTransformParameter( "TranslationZ" )
	print "translation:", tx, ty, tz
        error = (abs(tx - TRANSLATION_X) +
                 abs(ty - TRANSLATION_Y) +
                 abs(tz - TRANSLATION_Z))
        error = error/(TRANSLATION_X + TRANSLATION_Y + TRANSLATION_Z)
	print "error:", error
        correctResult = ( error < 0.5 )
        self.failUnless( correctResult, "Registration is incorrect")

        # registration.SetMetricType( VTK_METRIC_NORMALIZED_CROSS_CORRELATION )
        registration.SetMetricParameter( "AComponentExtent", 64 )
        registration.SetMetricParameter( "BComponentExtent", 64 )

        registration.SetTransformParameter( "CenterX", 31.5 )
        registration.SetTransformParameter( "CenterY", 31.5 )
        registration.SetTransformParameter( "CenterZ", 31.5 )
        registration.SetTransformParameter( "TranslationX", tx )
        registration.SetTransformParameter( "TranslationY", ty )
        registration.SetTransformParameter( "TranslationZ", tz )
        registration.SetTransformParameter( "RotationAngleX", 0.0 )
        registration.SetTransformParameter( "RotationAngleY", 0.0 )
        registration.SetTransformParameter( "RotationAngleZ", 0.0 )
        registration.SetTransformParameter( "IsotropicScale", 1.0 )
        registration.SetTransformParameter( "TranslationScale", 3.0 )
        registration.SetTransformParameter( "RotationScale", 0.0 )
        registration.SetTransformParameter( "ScalingScale", 0.0 )
 
        registration.UpdateRegistration()

        tx = registration.GetTransformParameter( "TranslationX" )
        ty = registration.GetTransformParameter( "TranslationY" )
        tz = registration.GetTransformParameter( "TranslationZ" )
        print "translation:", tx, ty, tz
        error = (abs(tx - TRANSLATION_X) +
                 abs(ty - TRANSLATION_Y) +
                 abs(tz - TRANSLATION_Z))
        error = error/(TRANSLATION_X + TRANSLATION_Y + TRANSLATION_Z)
        print "error:", error
        correctResult = ( error < 0.5 )
        self.failUnless( correctResult, "Registration is incorrect")

if __name__ == "__main__":
    Testing.main([(TestImageRegistration, 'test')])
