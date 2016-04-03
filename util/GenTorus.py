'''
Generates a torus based on parameters from the specified ini file
'''

import sys
import vtk
import numpy as np
from configobj import ConfigObj

def GenTorus(programArguments):
    ''' programArguments: ini file containing model parameters'''
    
    # Load relevant parameters from ini file
    conf = ConfigObj(programArguments)
    parameters = conf['Parameters']
    majorCirc = parameters['surfaceLength']
    minorCirc = parameters['surfaceWidth']
    thetaMesh = parameters['xMesh']

    # Minor and major radii
    r = float(minorCirc)/(2*np.pi)
    R = float(majorCirc)/(2*np.pi)

    # Mesh sizes
    thetaResolution = int(thetaMesh)    
    phiResolution = int(thetaResolution*(R/r))
    
    # Generate a torus 
    torusSource = vtk.vtkSuperquadricSource()
    torusSource.SetCenter(0.0, 0.0, 0.0)
    torusSource.SetScale(1.0, 1.0, 1.0)
    torusSource.SetToroidal(1) 
    torusSource.SetThetaRoundness(1)
    torusSource.SetPhiRoundness(1)

    # SuperquadricSource reverses phi and theta to be confusing - this is not an error
    torusSource.SetPhiResolution(thetaResolution)   
    torusSource.SetThetaResolution(phiResolution)

    # Don't change these!
    torusSource.SetSize(R + r)      
    torusSource.SetThickness(r/R)    

    # The quadric has nasty discontinuities from the way the edges are generated
    # so let's pass it though a CleanPolyDataFilter and merge any points which
    # are coincident, or very close. First convert quads into triangles
    tri = vtk.vtkTriangleFilter()
    tri.SetInputConnection(torusSource.GetOutputPort())
    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputConnection(tri.GetOutputPort())
    cleaner.SetTolerance(0.00005)
    cleaner.Update()
    
    outputFileName = "torus_R" + majorCirc + "_r" + minorCirc + "_mesh" + thetaMesh + ".vtp"

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(cleaner.GetOutput())
    writer.SetFileName(outputFileName)
    writer.Update()
    
    print "Saving output to file", outputFileName 

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print "Usage: " + sys.argv[0] + " <Program Arguments>"
    else:
        GenTorus(sys.argv[1])
