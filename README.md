CRD Model
=========

This project models reaction diffusion equations on different curved surfaces. 

Currently we have two models, the Fitz Hugh Nagumo neuron model and the Goldbeter smooth muscle cell model.

The Fitz Hugh Nagumo model contains an activator variable and an inhibitor variable. The activator variable diffuses over the surface, while the diffusion of the inhibitor variable is seen to be negligible and thus omitted. 
The Goldbeter model contains variables for the calcium concentration in the cytosol and concentration in the internal cell stores. Only the cytosolic calcium diffuses.

Both models are applied to the surface of a torus. A torus is chosen as it has areas of both negative and positive Gaussian curvature; this curvature affects the rate of diffusion. Hence the behaviour on a torus is different to the same equations modelled on a flat 2D surface.



How to Compile
--------------

The project depends on MPI and Sundials libraries.

Configure the project with CMake and specify the *out-of-source* build directory.

After the the project has been configured, it can be opened and compiled with
the target IDE, or in the case of Unix Makefile configuration simply run make in
the build directory. 

How to Run
----------

Modify the model parameters with the file FHNmodelArgs.ini or GoldbeterModelArgs.ini found in the data directory or supply your own ini file. This is included as a command line argument when running most scripts.


To run: 

    mpirun -np <<number of processes>> <<C program>> <<ini file>>

where C program is found in the src directory and depends on the surface and model needed. This outputs (Number of Variables)*(Number of processes) text files, where each text file corresponds to one variable in one subdomain. The rows correspond to each time step and the columns correspond to their grid point in the subdomain.


To plot as a 3D surface plot: 

    python <<plotting script>>  

where plotting script is plot_*.py found in util/FHNmodel or util/GoldbeterModel directory. This outputs a .png file for every timestep and combines them into .gif format.


To plot solutions on a torus, first generate a torus with the appropriate mesh size:

    python <<script>> <<ini file>>

where script is GenTorus.py, found in util directory. This outputs a .vtp file for use by MapOutputToTorus.py. This step can be omitted once the appropriate torus is generated.

Next map the output of FHNmodel_torus to the generated torus: 

    python <<script>> <<ini file>>
     
where script is MapOutputToTorus.py, found in util/FHNmodel or util/GoldbeterModel directory. This outputs multiple .vtp files called either FHN or GBstep_*.vtp where each file corresponds to a timestep, and an additional Paraview file called either FHNtimeSteps.pvd or GBtimeSteps.pvd, which contain the physiological timesteps of each step. This pvd file can be then opened in Paraview. It is recommended to use the Annotate Time Filter to see the time at each step.
