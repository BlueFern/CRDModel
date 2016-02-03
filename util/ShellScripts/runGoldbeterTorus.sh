#!/bin/bash

# Modify the paths as needed and make sure script is executable (chmod u+x)

# Run the C script using parameters in the ini file
time mpirun -np 4 ~/Documents/Research/CRDModel.x64/GoldbeterModel_torus ~/Documents/Research/CRDModel/data/GoldbeterModelArgs.ini

# Plot the solution on a 2D graph and combine into a video
time python ~/Documents/Research/CRDModel/util/GoldbeterModel/plot_GoldbeterModel_torus.py ~/Documents/Research/CRDModel/data/GoldbeterModelArgs.ini

# Map the solution onto a torus
#time python ~/Documents/Research/CRDModel/util/GoldbeterModel/MapOutputToTorus.py ~/Documents/Research/CRDModel/data/GoldbeterModelArgs.ini
