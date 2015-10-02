# Plot for periodic FHN model

# imports
import vtk
import sys
import numpy as np
import math

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
    # using the number of rows (ny), number of columns (nx)
    rc = (int(phi/(2*np.pi)*(ny-1)), int(theta/(2*np.pi)*(nx-1)) )

    return rc

def MapOutputToTorus(majorCirc, thetaMesh):
    ''' majorCirc = major circumference of the torus. Set to either 40 or 80
        thetaMesh = number of mesh points in theta direction '''

    # Minor radius of torus
    r = 20/(2*np.pi)
    
    # Major radius of torus
    R = float(majorCirc)/(2*np.pi)
    
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
       
       
    # load first processor's data, and determine total number of time steps
    data = np.loadtxt('FHNmodel_torus_u.000.txt', dtype=np.double)
    nt = np.shape(data)[0]
    # shape returns the dimensions of the array
    
    # create empty array for all solution data of the first variable
    results = np.zeros((nt,ny,nx))
    
    # insert first processor's data into results array
    istart = subdomains[0,0]
    iend = subdomains[0,1]
    jstart = subdomains[0,2]
    jend = subdomains[0,3]
    nxl = iend-istart+1
    nyl = jend-jstart+1
    for i in range(nt):
        results[i,jstart:jend+1,istart:iend+1] = np.reshape(data[i,:], (nyl,nxl))
        
    # iterate over remaining data files, inserting into output
    if (nprocs > 1):
        for isub in range(1,nprocs):
            data = np.loadtxt('FHNmodel_torus_u.' + repr(isub).zfill(3) + '.txt', dtype=np.double)
            # check that subdomain has correct number of time steps
            if (np.shape(data)[0] != nt):
                sys.exit('error: subdomain ' + isub + ' has an incorrect number of time steps')
            istart = subdomains[isub,0]
            iend = subdomains[isub,1]
            jstart = subdomains[isub,2]
            jend = subdomains[isub,3]
            nxl = iend-istart+1
            nyl = jend-jstart+1
            for i in range(nt):
                results[i,jstart:jend+1,istart:iend+1] = np.reshape(data[i,:], (nyl,nxl))
                
    for tstep in range(nt):
        
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
        
        activatorArray = vtk.vtkDoubleArray()
        activatorArray.SetName("Activator")
        
        # Iterate over all centres
        for cId in range(cellCentres.GetNumberOfPoints()):
            point = cellCentres.GetPoint(cId)
            
            rc = XYZtoRC(point,ny,nx,r,R)              
            
            result = results[tstep,rc[0],rc[1]]
            
            activatorArray.InsertNextValue(result)
    
        
        torus.GetCellData().SetScalars(activatorArray)
        
        outputFileName = "step_" + repr(tstep).zfill(3) + ".vtp"
        
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetInput(torus)
        writer.SetFileName(outputFileName)
        writer.Update()

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "Usage: " + sys.argv[0] + " <Major Circumference> <Theta Mesh Size>"
    else:
        MapOutputToTorus(sys.argv[1], sys.argv[2])
