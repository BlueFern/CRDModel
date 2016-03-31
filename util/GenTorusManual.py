'''
Generate a torus with manually specified R, r and mesh sizes
'''


import sys
import vtk
import numpy as np
from configobj import ConfigObj

# Specify the torus parameters
r = 5/(2*np.pi); # Minor radius
R = 6/(2*np.pi); # Major radius
thetaResolution = 200; # Mesh in theta direction
phiResolution = thetaResolution*(R/r); # Mesh in phi direction

torusSource = vtk.vtkSuperquadricSource()
torusSource.SetCenter(0.0, 0.0, 0.0)
torusSource.SetScale(1.0, 1.0, 1.0)
torusSource.SetToroidal(1)
torusSource.SetThetaRoundness(1)
torusSource.SetPhiRoundness(1)
torusSource.SetPhiResolution(thetaResolution)   # SuperquadricSource reverses phi and theta to be confusing - this is not an error

torusSource.SetSize(R+r)
torusSource.SetThickness(r)
torusSource.SetThetaResolution(phiResolution)

tri = vtk.vtkTriangleFilter()
tri.SetInputConnection(torusSource.GetOutputPort())

# The quadric has nasty discontinuities from the way the edges are generated
# so let's pass it though a CleanPolyDataFilter and merge any points which
# are coincident, or very close
cleaner = vtk.vtkCleanPolyData()
cleaner.SetInputConnection(tri.GetOutputPort())
cleaner.SetTolerance(0.00005)

cleaner.Update()

outputFileName = "torus_manual.vtp"
print "Saving output to file", outputFileName

writer = vtk.vtkXMLPolyDataWriter()
writer.SetInput(cleaner.GetOutput())
writer.SetFileName(outputFileName)
writer.Update()

