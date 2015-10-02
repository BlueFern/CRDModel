''' Plot the first variable u of the Fitz-Hugh Nagumo model on a torus, shown 
    as a flat surface with periodic boundaries '''

# imports
import sys
import numpy as np
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import os



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
        xmin = subd[6]
        xmax = subd[7]
    else:
        if ((subd[0] != nx) or (subd[1] != ny)):
            sys.exit("error: subdomain files incompatible (clean up and re-run test)")
    subdomains[i,:] = subd[2:6]

   
   
# First variable ****   
# load first processor's data, and determine total number of time steps
data = np.loadtxt('FHNmodel_torus_u.000.txt', dtype=np.double)
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

# determine extents of plots
maxtemp = 1.1*results.max()
mintemp = 0.9*results.min()


# generate plots of results
for tstep in range(nt):

    # set string constants for output plots, current time, mesh size
    pname = 'FHNmodel_torus_surf_u.' + repr(tstep).zfill(3) + '.png'
    tstr  = repr(tstep)
    nxstr = repr(nx)
    nystr = repr(ny)

    # set x and y meshgrid objects 
    xspan = np.linspace(xmin, xmax, nx)   
    yspan = np.linspace(0.0, 2.0*np.pi, ny)
    X,Y = np.meshgrid(xspan,yspan)

    # plot current solution as a surface, and save to disk
    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, results[tstep,:,:], rstride=1, cstride=1, 
                    cmap=cm.jet, vmin=mintemp, vmax=maxtemp, linewidth=0, antialiased=True, shade=True)
    
    ax.set_xlabel('theta')
    ax.set_ylabel('phi')

    ax.set_zlim((mintemp, maxtemp))
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    ax.view_init(90,90)
    #title('Torus (curved): u(x,y) at output ' + tstr + ', mesh = ' + nxstr + 'x' + nystr)
    title('Torus: u(theta, phi) at output ' + tstr + ', mesh = ' + nxstr + 'x' + nystr)
    savefig(pname)
    plt.close()


#print '\nConverting png files to animated gif (this maye take some time)...\n'
os.system("convert -delay 30 -loop 0 FHNmodel_torus_surf_u.*.png " + "FHNmodel_torus_outside_beta_09_R80.gif")
os.system("rm FHNmodel_torus_surf_u.*.png") # clean up files
os.system("animate FHNmodel_torus_outside_beta_09_R80.gif") # play animates gif






### Second variable
#
#
##load first processor's data, and determine total number of time steps
#data2 = np.loadtxt('FHNmodel_v.000.txt', dtype=np.double)
#nt = np.shape(data2)[0]
#
## create empty array for all solution data
#results = np.zeros((nt,ny,nx))
#
## insert first processor's data into results array
#istart = subdomains[0,0]
#iend = subdomains[0,1]
#jstart = subdomains[0,2]
#jend = subdomains[0,3]
#nxl = iend-istart+1
#nyl = jend-jstart+1
#for i in range(nt):
#    results[i,jstart:jend+1,istart:iend+1] = np.reshape(data2[i,:], (nyl,nxl))
#    
## iterate over remaining data files, inserting into output
#if (nprocs > 1):
#    for isub in range(1,nprocs):
#        data2 = np.loadtxt('FHNmodel_torus_v.' + repr(isub).zfill(3) + '.txt', dtype=np.double)
#        # check that subdomain has correct number of time steps
#        if (np.shape(data2)[0] != nt):
#            sys.exit('error: subdomain ' + isub + ' has an incorrect number of time steps')
#        istart = subdomains[isub,0]
#        iend = subdomains[isub,1]
#        jstart = subdomains[isub,2]
#        jend = subdomains[isub,3]
#        nxl = iend-istart+1
#        nyl = jend-jstart+1
#        for i in range(nt):
#            results[i,jstart:jend+1,istart:iend+1] = np.reshape(data2[i,:], (nyl,nxl))
#
## determine extents of plots
#maxtemp = 1.1*results.max()
#mintemp = 0.9*results.min()
#
## generate plots of results
#for tstep in range(nt):
#
#    # set string constants for output plots, current time, mesh size
#    pname = 'FHNmodel_surf_v.' + repr(tstep).zfill(3) + '.png'
#    tstr  = repr(tstep)
#    nxstr = repr(nx)
#    nystr = repr(ny)
#
#    # set x and y meshgrid objects
#    xspan = np.linspace(0.0, 2.0*np.pi, nx)
#    X,Y = np.meshgrid(xspan,yspan)
#
#    # plot current solution as a surface, and save to disk
#    fig = plt.figure(1)
#    ax = fig.add_subplot(111, projection='3d')
#    ax.plot_surface(X, Y, results[tstep,:,:], rstride=1, cstride=1, 
#                    cmap=cm.jet, vmin=mintemp, vmax=maxtemp, linewidth=0, antialiased=True, shade=True)
#    ax.set_xlabel('x')
#    ax.set_ylabel('y')
#    ax.set_zlim((mintemp, maxtemp))
#    plt.gca().invert_xaxis()
#    plt.gca().invert_yaxis()
#    ax.view_init(20,45)
#    title('v(x,y) at output ' + tstr + ', mesh = ' + nxstr + 'x' + nystr)
#    savefig(pname)
#    plt.close()





##### end of script #####
