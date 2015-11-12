# Show Gaussian curvature and coupling strength of the torus

# imports
import vtk
import sys
import numpy as np
import math
from configobj import ConfigObj

def XYZtoPT(xyz,r,R):
    '''Convert (x,y,z) to (phi,theta)'''
    
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    # Convert (x,y,z) coordinates to (phi, theta).
    phi = math.atan2(z,x) % (2*np.pi)
    
    if (np.sqrt(x*x + z*z) > R):    # i.e. if the point lies on the outer half of the torus
        theta = (np.arcsin(y/r)) % (2*np.pi)   
    else:
        theta = (np.pi -  np.arcsin(y/r)) % (2*np.pi)    
    
    return phi, theta
    
def PTtoETA(phi, theta, r, R):
    '''Convert (phi, theta) to eta, theta_i, a - alternate toroidal coordinates
        used in the Coupling Strength equation'''

    a = np.sqrt(math.pow(R,2) - math.pow(r,2))    
    
    eta = math.atanh(a/R)
    
    if theta >= 0 and theta <= np.pi: 
        theta_i = math.acos( R/r - math.pow(a,2)/(r*(R+r*np.cos(theta))))  
    else:
        theta_i = -math.acos( R/r - math.pow(a,2)/(r*(R+r*np.cos(theta))))    
    
    return eta, theta_i, a
    

def GenCurvatureCoupling(programArguments):
    ''' programArguments: ini file containing model parameters'''
    
    # Load relevant parameters from ini file
    conf = ConfigObj(programArguments)
    parameters = conf['Parameters']
    majorCirc = parameters['majorCirc']
    thetaMesh = parameters['thetaMesh']

    # Minor radius of torus
    r = 20/(2*np.pi)
    
    # Major radius of torus
    R = float(majorCirc)/(2*np.pi)
    
    
    # Read geometry from disk
    torusReader = vtk.vtkXMLPolyDataReader()
    torusReader.SetFileName("torus_R" + majorCirc + "_mesh" + thetaMesh + ".vtp")
    torusReader.Update()
    
    torus = torusReader.GetOutput()
    
    # Obtain cell centres
    cellCentresFilter = vtk.vtkCellCenters()
    cellCentresFilter.SetInput(torus)
    cellCentresFilter.Update()

    cellCentres = cellCentresFilter.GetOutput()
    
    gaussianCurvature = vtk.vtkDoubleArray()
    gaussianCurvature.SetName("Gaussian Curvature")
    
    couplingStrength = vtk.vtkDoubleArray()
    couplingStrength.SetName("Coupling Strength")
    
    # Iterate over all centres
    for cId in range(cellCentres.GetNumberOfPoints()):
        point = cellCentres.GetPoint(cId)
        
        phi, theta = XYZtoPT(point,r,R)              
        
        # Gaussian curvature
        resultG = np.cos(theta)/(r*(R+r*np.cos(theta)))
        
        #Convert phi, theta        
        eta, theta_i, a = PTtoETA(phi, theta, r, R)
        
        # Coupling strength
        resultC = 10*math.pow((math.cosh(eta) - np.cos(theta_i)),2)/math.pow(a,2)   # Not sure why 10 is there but is necessary
        
        gaussianCurvature.InsertNextValue(resultG)
        couplingStrength.InsertNextValue(resultC)
    
    torus.GetCellData().SetScalars(gaussianCurvature)
    torus.GetCellData().AddArray(couplingStrength)  
    
    outputFileName = "CurvatureCoupling_torus_R" + majorCirc + "_mesh" + thetaMesh + ".vtp"
    print "Saving output to file", outputFileName 
    
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInput(torus)
    writer.SetFileName(outputFileName)
    writer.Update() 
    
    
    
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print "Usage: " + sys.argv[0] + " <Program Arguments>"
    else:
        GenCurvatureCoupling(sys.argv[1])