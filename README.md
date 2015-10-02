CRD Model
=========

This project models reaction diffusion equations on different curved surfaces. 

Currently the FitzHugh Nagumo neuron model is applied to the surface of a torus. This model contains an activator variable and an inhibitor variable. The activator variable diffuses over the surface, while the diffusion of the inhibitor variable is seen to be negligible and thus omitted. 

A torus is chosen as it has areas of both negative and positive Gaussian curvature; this curvature affects the rate of diffusion of the activator variable. Hence the behaviour on a torus is different to the same equations modelled on a flat 2D surface.

How to Compile
--------------

The project depends on MPI and Sundials libraries.

Configure the project with CMake and specify the *out-of-source* build directory.

After the the project has been configured, it can be opened and compiled with
the target IDE, or in the case of Unix Makefile configuration simply run make in
the build directory. 

How to Run
----------

To do.....

