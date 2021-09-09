import numpy as np

"""Input parameters"""
A_b = 60  # [m^2] Air cushion area
h_b = 0.4  # [m] Distance from waterline to the roof of the air cushion
M = 25e3  # [kg] Total mass of the vessel
p_0 = 3500  # [Pa] Excess air pressure in air cushion
p_a = 101325  # [Pa] Atmospheric pressure
gamma = 1.4  # [-] Ratio of specific heat for air
Q_0 = 50  # [m^3/s] Mean volumetric flow of air fan
dQdp_0 = -7.782e-3  # [(m^3s^-1)/(Pa)] Slope of fan characteristic curve estimated from Figure 5.6 in Faltinsen High-Speed

'''
Equivalent mass, damping and spring coefficients
m*ddx + c*dx + k*x = F(t)
'''
m = M*h_b/gamma/(p_0 + p_a)  # [m^2 s^2] Mass equivalent coefficient
c = M/A_b/p_0 * (0.5*Q_0 - dQdp_0*p_0)  # [m^2 s] Damping equivalent coefficient
k = A_b  # [m^2] Stiffness equivalent coefficient

'''Calculating natural frequency and damping ratio'''

omega_0 = np.sqrt(k/m)  # [rad/s] Natural frequency of uniform pressure oscillations
c_cr = 2 * np.sqrt(k*m)  # [m^2 s] Critical damping of the system
zeta = c/c_cr  # [-] Damping ratio of the system

'''Print results'''

print('Natural frequency:\t', round(omega_0, 2), '\t[rad/s]')
print('Natural frequency:\t', round(omega_0/2*np.pi, 2), '\t[Hz]')
print('Damping ratio:\t\t', round(zeta, 3), '\t[-]')

