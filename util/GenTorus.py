import sys
import vtk
import numpy as np
from configobj import ConfigObj

def GenTorus(programArguments):
    ''' programArguments: ini file containing model parameters'''
    
    # Load relevant parameters from ini file
    conf = ConfigObj(programArguments)
    parameters = conf['Parameters']
    majorCirc = parameters['majorCirc']
    thetaMesh = parameters['thetaMesh']
    
    thetaResolution = int(thetaMesh)    
    
    torusSource = vtk.vtkSuperquadricSource()
    torusSource.SetCenter(0.0, 0.0, 0.0)
    torusSource.SetScale(1.0, 1.0, 1.0)
    torusSource.SetToroidal(1) 
    torusSource.SetThetaRoundness (1)
    torusSource.SetPhiRoundness (1)
    torusSource.SetPhiResolution(thetaResolution)   # SuperquadricSource reverses phi and theta to be confusing - this is not an error

    torusSource.SetSize(float(majorCirc)/(2*np.pi)+20/(2*np.pi))      
    torusSource.SetThickness (20/float(majorCirc))    
    torusSource.SetThetaResolution(thetaResolution*(float(majorCirc)/20))
    #torusSource.SetThetaResolution(thetaResolution)
    
    tri = vtk.vtkTriangleFilter()
    tri.SetInputConnection(torusSource.GetOutputPort())
     
    # The quadric has nasty discontinuities from the way the edges are generated
    # so let's pass it though a CleanPolyDataFilter and merge any points which
    # are coincident, or very close
    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputConnection(tri.GetOutputPort())
    cleaner.SetTolerance(0.00005)
    
    cleaner.Update()
    
    #outputFileName = "torus_R" + majorCirc + "_tmesh" + thetaMesh + "_pmesh" + thetaMesh + ".vtp"
    outputFileName = "torus_R" + majorCirc + "_mesh" + thetaMesh + ".vtp"
    print "Saving output to file", outputFileName 
    
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInput(cleaner.GetOutput())
    writer.SetFileName(outputFileName)
    writer.Update()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print "Usage: " + sys.argv[0] + " <Program Arguments>"
    else:
        GenTorus(sys.argv[1])
