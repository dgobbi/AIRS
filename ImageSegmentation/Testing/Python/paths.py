"""
Importing this module will add to your module search path
"""

import sys, os.path

## # paths.prefix is the path to the atamai directory
## mod = sys.modules[__name__]
## file = os.path.abspath(mod.__file__)
## exampledir = os.path.split(file)[0]
## prefix = os.path.split(exampledir)[0]
## del mod
## del file
## del exampledir

path,name = os.path.split(__file__)

sys.path.insert(0,os.path.join(path,"atamai/classes"))
sys.path.insert(0,os.path.join(path,"atamai/modules"))
sys.path.insert(0,os.path.join(path,"atamai/minc"))
sys.path.insert(0,os.path.join(path,"atamai/misc"))
sys.path.insert(0,os.path.join(path,"atamai/dicom"))
##sys.path.insert(0,os.path.join(path,"vtkMultiIO"))

