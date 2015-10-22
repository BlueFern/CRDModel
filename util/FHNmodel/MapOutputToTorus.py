# Plot for periodic FHN model

# imports
import vtk
import os
import sys
import lxml
import lxml.etree
import numpy as np
import argparse
import math
from configobj import ConfigObj

def XYZtoRC(xyz,ny,nx,r,R):
    '''Convert (x,y,z) to (phi,theta) to (row,column) of results array'''
    
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    # Convert (x,y,z) coordinates to (phi, theta).
    phi = math.atan2(z,x) % (2*np.pi)
    
    if (np.sqrt(x*x + z*z) > R):    # i.e. if the point lies on the outer half of the torus
        theta = (np.arcsin(y/r)) % (2*np.pi)   
    else:
        theta = (np.pi -  np.arcsin(y/r)) % (2*np.pi)

    # Finds the row and column of the results array that match the given coords
    # using the number of rows (ny) and number of columns (nx)
    rc = (int(phi/(2*np.pi)*(ny-1)), int(theta/(2*np.pi)*(nx-1)) )

    return phi, theta, rc
    

if __name__ == '__main__':

    # Obtain ini file
    parser = argparse.ArgumentParser()
    parser.add_argument("configFile", help="config file path")
    args = parser.parse_args()

    # Load relevant parameters from ini file
    conf = ConfigObj(args.configFile)
    parameters = conf['Parameters']
    systemParameters = conf['System']
    majorCirc = parameters['majorCirc']
    thetaMesh = parameters['thetaMesh']
    tFinal = parameters['tFinal']
    includeAllVars = int(systemParameters['includeAllVars'])
    varyBeta = int(systemParameters['varyBeta'])

    # Minor radius of torus
    r = 20/(2*np.pi)
    
    # Major radius of torus
    R = float(majorCirc)/(2*np.pi)

    BETAMIN = 0.7
    BETAMAX = 1.7

    # Hopf bifurcation position
    if (varyBeta == 1): Hopf = (1.0 - BETAMIN)*2*np.pi/(BETAMAX - BETAMIN)
    
    # determine the number of MPI processes used
    nprocs=1
    for i in range(1000):
        sname = 'FHNmodel_torus_subdomain.' + repr(i).zfill(3) + '.txt'
        try:
            f = open(sname,'r')
            f.close()
        except IOError:
            nprocs = i
            break
    
    # load subdomain information, store in table
    subdomains = np.zeros((nprocs,4), dtype=np.int)
    for i in range(nprocs):
        sname = 'FHNmodel_torus_subdomain.' + repr(i).zfill(3) + '.txt'
        subd = np.loadtxt(sname, dtype=np.float)
        if (i == 0):
            nx = int(subd[0])
            ny = int(subd[1])
        else:
            if ((subd[0] != nx) or (subd[1] != ny)):
                sys.exit("error: subdomain files incompatible (clean up and re-run test)")
        subdomains[i,:] = subd[2:6]
       
       
    # load first processor's data for variable U and V
    dataU = np.loadtxt('FHNmodel_torus_u.000.txt', dtype=np.double)
    if (includeAllVars == 1):  dataV = np.loadtxt('FHNmodel_torus_v.000.txt', dtype=np.double)
    
    # determine total number of time steps
    nt = np.shape(dataU)[0]
    
    # create empty array for all solution data of the variables U and V
    resultsU = np.zeros((nt,ny,nx))
    if (includeAllVars == 1): resultsV = np.zeros((nt,ny,nx))
    
    # insert first processor's data into results array
    istart = subdomains[0,0]
    iend = subdomains[0,1]
    jstart = subdomains[0,2]
    jend = subdomains[0,3]
    nxl = iend-istart+1
    nyl = jend-jstart+1
    
    for i in range(nt):
        resultsU[i,jstart:jend+1,istart:iend+1] = np.reshape(dataU[i,:], (nyl,nxl))
        if (includeAllVars == 1): resultsV[i,jstart:jend+1,istart:iend+1] = np.reshape(dataU[i,:], (nyl,nxl))
        
    # iterate over remaining data files, inserting into output
    if (nprocs > 1):
        
        for isub in range(1,nprocs):
            dataU = np.loadtxt('FHNmodel_torus_u.' + repr(isub).zfill(3) + '.txt', dtype=np.double)
            if (includeAllVars == 1): dataV = np.loadtxt('FHNmodel_torus_v.' + repr(isub).zfill(3) + '.txt', dtype=np.double)
            
            # check that subdomain has correct number of time steps
            if (np.shape(dataU)[0] != nt):
                sys.exit('error: subdomain ' + isub + ' has an incorrect number of time steps')
                
            istart = subdomains[isub,0]
            iend = subdomains[isub,1]
            jstart = subdomains[isub,2]
            jend = subdomains[isub,3]
            nxl = iend-istart+1
            nyl = jend-jstart+1
            
            for i in range(nt):
                resultsU[i,jstart:jend+1,istart:iend+1] = np.reshape(dataU[i,:], (nyl,nxl))
                if (includeAllVars == 1): resultsV[i,jstart:jend+1,istart:iend+1] = np.reshape(dataV[i,:], (nyl,nxl))

    tStepDict = dict()

    # Create directory for each time step output vtp file
    os.system("mkdir FHNstep")

    for tstep in range(nt):
        
        # Get the time at tstep
        time = ((float(tstep)/float(nt)))*float(tFinal) # get time of output
        
        # Read geometry from disk
        torusReader = vtk.vtkXMLPolyDataReader()
        #torusReader.SetFileName("torus_R" + majorCirc + "_tmesh" + thetaMesh + "_pmesh" + thetaMesh + ".vtp")
        torusReader.SetFileName("torus_R" + majorCirc + "_mesh" + thetaMesh + ".vtp")
        torusReader.Update()
        
        torus = torusReader.GetOutput()
        
        # Obtain cell centres
        cellCentresFilter = vtk.vtkCellCenters()
        cellCentresFilter.SetInput(torus)
        cellCentresFilter.Update()
    
        cellCentres = cellCentresFilter.GetOutput()
        
        activatorArray = vtk.vtkDoubleArray()
        activatorArray.SetName("Activator")

        if (varyBeta == 1):
            HopfArray = vtk.vtkDoubleArray()
            HopfArray.SetName("Hopf Bifurcations")
        
        if (includeAllVars == 1):
            inhibitorArray = vtk.vtkDoubleArray()
            inhibitorArray.SetName("Inhibitor")
        
        # Iterate over all centres
        for cId in range(cellCentres.GetNumberOfPoints()):
            point = cellCentres.GetPoint(cId)
            
            phi, theta, rc = XYZtoRC(point,ny,nx,r,R)
            
            resultU = resultsU[tstep,rc[0],rc[1]]
            if (includeAllVars == 1): resultV = resultsV[tstep,rc[0],rc[1]]

            # Set to 1 if it's at one of the Hopf bifurcation, otherwise 0
            if (varyBeta == 1):
                if abs(phi - Hopf) < 0.01:
                    resultHopf = 1
                else:
                    resultHopf = 0

            activatorArray.InsertNextValue(resultU)
            if (varyBeta == 1): HopfArray.InsertNextValue(resultHopf)
            if (includeAllVars == 1): inhibitorArray.InsertNextValue(resultV)
        
        torus.GetCellData().SetScalars(activatorArray)  # default variable
        if (varyBeta == 1): torus.GetCellData().AddArray(HopfArray) # Hopf bifurcations
        if (includeAllVars == 1): torus.GetCellData().AddArray(inhibitorArray)    # other variable(s)


        outputFileName = "FHNstep/FHNstep_" + repr(tstep).zfill(3) + ".vtp"

        tStepDict[time] = outputFileName
        
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetInput(torus)
        writer.SetFileName(outputFileName)
        writer.Update()

    # Create new root of the XML tree.
    vtkFileElem = lxml.etree.Element('VTKFile', type="Collection", version='0.1', byte_order='LittleEndian', compressor='vtkZLibDataCompressor')

    # Create new XML document.
    pvdDoc = lxml.etree.ElementTree(vtkFileElem)

    # Add 'Collections' element, as a subelement of the 'VTKFile'.
    collectionElem = lxml.etree.SubElement(vtkFileElem, 'Collection')

    for time in tStepDict:
        tstr = repr(float("{0:.1f}".format(time)))  # convert to string with 1 decimal place
        # Add 'DataSet' element as a subelement of ... whatever.
        dataElem = lxml.etree.SubElement(collectionElem, 'DataSet', timestep=tstr, group='', part='0', file=tStepDict[time])

    # Print the tree to the console.
    outXMLFile = open('FHNtimeSteps.pvd', 'w')
    pvdDoc.write(outXMLFile, encoding='iso-8859-1', xml_declaration=True, pretty_print=True, method = 'xml')


