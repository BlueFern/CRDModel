#!/bin/bash

# Modify the following paths and .bashrc file as needed and make sure script is executable (chmod u+x)

# Run the C script using parameters in the ini file
time mpirun -np 4 ~/Documents/Research/CRDModel.x64/FHNmodel_torus ~/Documents/Research/CRDModel/data/FHNmodelArgs.ini

# Plot the solution on a 2D graph and combine into a video
time python ~/Documents/Research/CRDModel/util/FHNmodel/plot_FHNmodel_torus.py ~/Documents/Research/CRDModel/data/FHNmodelArgs.ini

# Generate a torus
time python ~/Documents/Research/CRDModel/util/GenTorus.py ~/Documents/Research/CRDModel/data/FHNmodelArgs.ini 

# Map the solution onto a torus
time python ~/Documents/Research/CRDModel/util/FHNmodel/MapOutputToTorus.py ~/Documents/Research/CRDModel/data/FHNmodelArgs.ini
