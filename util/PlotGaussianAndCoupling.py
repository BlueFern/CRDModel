'''
Plots the Gaussian curvature G(theta) and the coupling strength C(theta) against
theta, for a flat surface, strongly curved and weakly curved torus.
'''


import matplotlib.pyplot as plt
import numpy as np
import math

# Gaussian Curvature
def G(t, r, R):
    return np.cos(t) / (r * (R + r * np.cos(t)))

# Coupling Strength - see Kneer et al. (2014) or thesis for details
def C(t, r, R):
    a = np.sqrt(R**2 - r**2)
    eta = math.atanh(a/R)
    theta2 = np.arccos(R/r - a**2 / (r * (R + r * np.cos(t))))
    return 10*(math.cosh(eta) - np.cos(theta2))**2 / a**2

r = 20 / (2* np.pi)
R1 = 80 / (2* np.pi)
R2 = 40 / (2* np.pi)
t = np.linspace(0 , 2*np.pi, 500)

# Plot Gaussian curvature
plt.xlabel('$\\theta$',fontsize=20)
plt.title('Gaussian Curvature $G(\\theta)$',fontsize=20)
plt.xlim(0, 2*np.pi)
plt.plot(t, G(t,r,R1), label='$R = 80 / 2 \pi$')
plt.plot(t, G(t,r,R2), "--", label='$R = 40 / 2 \pi$')
plt.plot((0,2*np.pi), (0,0), ':',label='Flat')
plt.legend(loc=4,fontsize=20)
plt.show()

# Plot Coupling strength
# plt.xlabel('$\\theta$',fontsize=20)
# plt.title('Coupling Strength $C(\\theta)$',fontsize=20)
# plt.xlim(0, 2*np.pi)
# plt.plot(t, C(t,r,R1), label='$R = 80 / 2 \pi$')
# plt.plot(t, C(t,r,R2), "--", label='$R = 40 / 2 \pi$')
# plt.plot((0,2*np.pi), (1,1), ':',label='Flat')
# plt.legend(loc=1,fontsize=20)
# plt.show()
