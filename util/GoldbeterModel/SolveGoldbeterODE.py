#!/usr/bin/env python

import numpy as np
from scipy import integrate
from matplotlib.pylab import *
import argparse
from configobj import ConfigObj

# Declare global beta parameter, obtained later from ini file
beta = 0.0

def goldbeter(t, y):

    v_0 = 1
    k = 10
    k_f = 1
    v_1 = 7.3
    V_M2 = 65
    V_M3 = 500
    K_2 = 1
    K_R = 2
    K_A = 0.9
    m = 2
    n = 2
    p = 4

    # Assign some variables for convenience of notation
    Z = y[0]
    Y = y[1]

    # Algebraic equations
    v_2 = V_M2 * Z**n / ( K_2**n + Z**n )
    v_3 = V_M3 * Y**m * Z**p / ( (K_R**m + Y**m) * (K_A**p + Z**p) )

    # Output from ODE function must be a COLUMN vector, with n rows
    n = len(y)      # 2: implies we have two ODEs
    dydt = np.zeros((n,1))
    dydt[0] = v_0 + v_1*beta - v_2 + v_3 + k_f*Y - k*Z
    dydt[1] = v_2 - v_3 - k_f*Y
    return dydt

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("beta", type=float, help="beta parameter")
    args = parser.parse_args()

    beta = args.beta

    # Start by specifying the integrator:
    # use ``vode`` with "backward differentiation formula"
    r = integrate.ode(goldbeter).set_integrator('vode', method='bdf')

    # Set the time range: go backwards in time to unstable FP if beta is in the oscillatory region
    if beta < 0.774 and beta > 0.289:
        t_start = 0.0
        t_final = -50.0
        delta_t = -0.1
    else:
        t_start = 0.0
        t_final = 50.0
        delta_t = 0.1

    # Number of time steps: 1 extra for initial condition
    num_steps = np.floor((t_final - t_start)/delta_t) + 1

    # Set initial condition(s): for integrating variable and time!
    Z_t_zero = 0.6
    Y_t_zero = 1.2
    r.set_initial_value([Z_t_zero, Y_t_zero], t_start)

    # Additional Python step: create vectors to store trajectories
    t = np.zeros((num_steps, 1))
    Z = np.zeros((num_steps, 1))
    Y = np.zeros((num_steps, 1))
    t[0] = t_start
    Z[0] = Z_t_zero
    Y[0] = Y_t_zero

    # Integrate the ODE(s) across each delta_t timestep
    k = 1
    while r.successful() and k < num_steps:
        r.integrate(r.t + delta_t)

        # Store the results to plot later
        t[k] = r.t
        Z[k] = r.y[0]
        Y[k] = r.y[1]
        k += 1

    # All done!  Plot the trajectories in two separate plots:
    fig = figure()
    ax1 = subplot(211)
    ax1.plot(t, Z)
    ax1.set_xlim(t_start, t_final)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Z')

    ax2 = plt.subplot(212)
    ax2.plot(t, Y)
    ax2.set_xlim(t_start, t_final)
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Y')

    print Z[-1], Y[-1]