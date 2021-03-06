''' Plot the first variable U of the FitzHugh Nagumo model on a flat surface
    with periodic boundaries '''

# imports
import sys
import numpy as np
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import os
from configobj import ConfigObj

def plot_FHNmodel_flat(programArguments):
    
    # Load relevant parameters from ini file
    conf = ConfigObj(programArguments)
    parameters = conf['Parameters']
    systemParameters = conf['System']
    surfaceLength = parameters['surfaceLength']
    surfaceWidth = parameters['surfaceWidth']
    beta = parameters['beta']
    tFinal = parameters['tFinal']
    betaMin = float(parameters['betaMin'])
    betaMax = float(parameters['betaMax'])
    varyBeta = int(systemParameters['varyBeta'])
    
    # determine the number of MPI processes used
    nprocs=1
    for i in range(1000):
        sname = 'FHNmodel_flat_subdomain.' + repr(i).zfill(3) + '.txt'
        try:
            f = open(sname,'r')
            f.close()
        except IOError:
            nprocs = i
            break
        
    # load subdomain information, store in table
    subdomains = np.zeros((nprocs,4), dtype=np.int)
    for i in range(nprocs):
        sname = 'FHNmodel_flat_subdomain.' + repr(i).zfill(3) + '.txt'
        subd = np.loadtxt(sname, dtype=np.float)
        if i == 0:
            nx = int(subd[0])
            ny = int(subd[1])
        else:
            if ((subd[0] != nx) or (subd[1] != ny)):
                sys.exit("error: subdomain files incompatible (clean up and re-run test)")
        subdomains[i,:] = subd[2:6]
       
       
    # First variable ****   
    # load first processor's data, and determine total number of time steps
    data = np.loadtxt('FHNmodel_flat_u.000.txt', dtype=np.double)
    nt = np.shape(data)[0]
    # shape returns the dimensions of the array
    
    # create empty array for all solution data
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
            data = np.loadtxt('FHNmodel_flat_u.' + repr(isub).zfill(3) + '.txt', dtype=np.double)
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
    
    # determine extents of plots
    maxtemp = results.max()
    mintemp = results.min()

    if varyBeta == 1:
        # Obtain location of Hopf (at beta = 1) by inverse of beta = BETAMIN + (BETAMAX - BETAMIN)/(surfacelength)*phi
        Hopf = (1.0 - betaMin)*float(surfaceLength)/(betaMax - betaMin)

    # Create subdirectory for the png files
    os.system("mkdir -p png")

    # generate plots of results
    for tstep in range(nt):
    
        # set string constants for output plots, current time, mesh size
        if varyBeta == 0:
            pname = 'png/FHNmodel_flat_Z.beta' + beta + '.' + repr(tstep).zfill(3) + '.png'
        else:
            pname = 'png/FHNmodel_flat_Z.varyBeta_linear' + repr(tstep).zfill(3) + '.png'

        time = ((float(tstep)/float(nt)))*float(tFinal) # get time of output
        tstr = repr(float("{0:.1f}".format(time)))  # convert to string with 1 decimal place
        nxstr = repr(nx)
        nystr = repr(ny)
    
        # set x and y meshgrid objects
        xspan = np.linspace(0.0, float(surfaceWidth), nx)
        yspan = np.linspace(0.0, float(surfaceLength), ny)
        X,Y = np.meshgrid(xspan,yspan)
    
        # plot current solution as a surface, and save to disk
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        img = ax.imshow(results[tstep,:,:], extent=[X.min(), X.max(), Y.min(), Y.max()], cmap='jet', aspect='auto', vmin=mintemp, vmax=maxtemp, origin='lower')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        fig.colorbar(img)

        if varyBeta == 1:
            plt.axhline(y=Hopf, color = 'r', linewidth=1, linestyle='dashed')

        time = ((float(tstep)/float(nt)))*float(tFinal) # get time of output
        tstr = repr(float("{0:.1f}".format(time)))  # convert to string with 1 decimal place
        nxstr = repr(nx)
        nystr = repr(ny)
        title('Flat: u(x,y) at t = ' + tstr + ', mesh = ' + nxstr + 'x' + nystr)
        savefig(pname, dpi=150)
        plt.close()
    
    
    # Convert png to video mp4
    if varyBeta == 0:
        os.system("ffmpeg -r 6 -i png/FHNmodel_flat_Z.beta" + beta + ".%03d.png FHNmodel_flat_Z.beta" + beta + ".mp4")
      # os.system("rm GoldbeterModel_flat_Z.*.png") # clean up files
    else:
        os.system("ffmpeg -r 6 -i png/FHNmodel_flat_Z.varyBeta_linear%03d.png FHNmodel_flat_Z.varyBeta_linear.mp4")
      # os.system("rm GoldbeterModel_flat_Z.varyBeta*.png") # clean up files

    print "\npng to mp4 conversion complete!"

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print "Usage: " + sys.argv[0] + " <Program Arguments>"
    else:
        plot_FHNmodel_flat(sys.argv[1])

