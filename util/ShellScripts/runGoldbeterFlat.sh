#!/bin/bash

# Modify the following paths and .bashrc file as needed and make sure script is executable (chmod u+x)

# Run the C script using parameters in the ini file
time mpirun -np 4 ~/Documents/Research/CRDModel.x64/GoldbeterModel_flat ~/Documents/Research/CRDModel/data/GoldbeterModelArgs.ini

# Plot the solution on a 2D graph and combine into a video
time python ~/Documents/Research/CRDModel/util/GoldbeterModel/plot_GoldbeterModel_flat.py ~/Documents/Research/CRDModel/data/GoldbeterModelArgs.ini
