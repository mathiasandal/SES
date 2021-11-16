from veres import read_veres_input, iterate_natural_frequencies
from air_cushion import read_fan_characteristics, air_cushion_area, interpolate_fan_characteristics, damping_matrix_air_cushion, stiffness_matrix_air_cushion
from mass_matrix import create_mass_matrix
from Wave_response_utilities import solve_eq_motion_steady_state, decouple_matrix, add_row_and_column

import matplotlib.pyplot as plt
import scipy.linalg as la
import numpy as np


# Read input from .re7 an .re8 files
path_veres = 'Input files//Veres input files//22kn//1-15s periods'
A_h, B_h, C_h, F_ex_real, F_ex_im, VEL, HEAD, FREQ, XMTN, ZMTN = read_veres_input(path_veres)

# Read fan characteristics
Q, P, rpm = read_fan_characteristics('Input files//fan characteristics//fan characteristics.csv', '1800rpm')

# Air cushion input variables
l_rect = 12  # [m] length of the rectangular part of the air cushion
l_tri = 6  # [m] length of the triangular part of the air cushion
b_c = 3.4  # [m] beam of the air cushion

h = 0.5  # [m] mean height between waterline and hull inside air cushion
z_c = -0.5 * h  # [m] vertical centroid of the air cushion relative to the ShipX coordinate system

p_0 = 3500  # [Pa] excess pressure in air cushion at equilibrium

S_0c, x_c = air_cushion_area(l_rect, l_tri, b_c)

Q_0, dQdp_0 = interpolate_fan_characteristics(p_0, P, Q)


# Create damping and stiffness matrix from air cushion #TODO: input should be difference between VERES coordinate system and cushion centroid
B_c = damping_matrix_air_cushion(S_0c, x_c - XMTN, h, p_0)  # Damping
C_c = stiffness_matrix_air_cushion(S_0c, h, x_c - XMTN, z_c - ZMTN, Q_0, dQdp_0, p_0)  # Stiffness

# Create mass matrix
# main dimensions of BBGreen
beam = 6  # [m] beam of BBGreen
Lpp = 19.4  # [m] L_pp of BBGreen
total_mass = 25.6e3  # [kg] total mass of the vessel
r44 = 0.35 * beam  # [m] radii of gyration in roll
r55 = 0.25 * Lpp  # [m] radii of gyration in pitch
r66 = 0.27 * Lpp  # [m] radii of gyration in yaw

# Creates mass matrix
M = create_mass_matrix(total_mass, r44, r55, r66)

# Collect stiffness matrices
C = add_row_and_column(C_h[0, 0, 0, :, :]) + C_c

# Decouple matrices to (3x3) containing heave, pitch and cushion pressure
decouple = False
if decouple:
    M = decouple_matrix(M, [2, 4, 6])
    C = decouple_matrix(C, [2, 4, 6])

# Compute natural frequencies and eigenmodes
# Iterate through frequencies to find the true natural frequencies

nat_frequencies_squared, eigen_modes, encounter_frequencies = iterate_natural_frequencies(FREQ, VEL[0], HEAD[0], A_h, M, C)

nat_frequencies = np.power(abs(nat_frequencies_squared), 0.5)

print(nat_frequencies)
eigen_modes = np.array(eigen_modes)
print(eigen_modes)
