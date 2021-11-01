from veres import read_veres_input
from air_cushion import read_fan_characteristics, air_cushion_area, interpolate_fan_characteristics, damping_matrix_air_cushion, stiffness_matrix_air_cushion
from mass_matrix import create_mass_matrix
from Wave_response_utilities import solve_eq_motion_steady_state, decouple_matrix

import matplotlib.pyplot as plt
import scipy.linalg as la
import numpy as np


# Read input from .re7 an .re8 files
path_veres = 'Input files//Veres input files'
A_h, B_h, C_h, F_ex_real, F_ex_im, VEL, HEAD, FREQ = read_veres_input(path_veres)

# Read fan characteristics
Q, P, rpm = read_fan_characteristics('Input files//fan characteristics//fan characteristics.csv', '1800rpm')

# Air cushion input variables
l_rect = 9  # [m] length of the rectangular part of the air cushion
l_tri = 10  # [m] length of the triangular part of the air cushion
b_c = 4  # [m] beam of the air cushion

h = 0.4  # [m] mean height between waterline and hull inside air cushion
z_c = 0.5 * h  # [m] vertical centroid of the air cushion relative to the ShipX coordinate system

p_0 = 3500  # [Pa] excess pressure in air cushion at equilibrium

S_0c, x_c = air_cushion_area(l_rect, l_tri, b_c)

Q_0, dQdp_0 = interpolate_fan_characteristics(p_0, P, Q)

# Create damping and stiffness matrix from air cushion
B_c = damping_matrix_air_cushion(S_0c, x_c, h, p_0)  # Damping
C_c = stiffness_matrix_air_cushion(S_0c, h, x_c, z_c, Q_0, dQdp_0, p_0)  # Stiffness

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
C = C_h + C_c

# Decouple matrices to (3x3) containing heave, pitch and cushion pressure


# Compute natural frequencies and eigenmodes
# Iterate through frequencies to find the true natural frequencies

