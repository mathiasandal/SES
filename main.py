from veres import read_veres_input, iterate_natural_frequencies, compute_RAOs
from air_cushion import read_fan_characteristics, air_cushion_area, interpolate_fan_characteristics, damping_matrix_air_cushion, stiffness_matrix_air_cushion
from mass_matrix import create_mass_matrix
from Wave_response_utilities import solve_eq_motion_steady_state, decouple_matrix, add_row_and_column

import matplotlib.pyplot as plt
import scipy.linalg as la
import numpy as np


# Read input from .re7 an .re8 files
path_veres = 'Input files//Veres input files//22kn//1-15s periods'
A_h, B_h, C_h, F_ex_real, F_ex_imag, VEL, HEAD, FREQ, XMTN, ZMTN = read_veres_input(path_veres)

# Read fan characteristics
Q, P, rpm = read_fan_characteristics('Input files//fan characteristics//fan characteristics.csv', '1800rpm')

# Air cushion input variables
l_rect = 12  # [m] length of the rectangular part of the air cushion
l_tri = 6  # [m] length of the triangular part of the air cushion
b_c = 3.4  # [m] beam of the air cushion

h = 0.64  # [m] mean height between waterline(baseline) and hull inside air cushion at AP
z_c = 0.5 * h  # [m] vertical centroid of the air cushion relative to the ShipX coordinate system

p_0 = 3500  # [Pa] excess pressure in air cushion at equilibrium

S_0c, x_c = air_cushion_area(l_rect, l_tri, b_c)  # Computes air cushion area properties

Q_0, dQdp_0 = interpolate_fan_characteristics(p_0, P, Q)  # Interpolates fan characteristics


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
nat_frequencies_squared, eigen_modes, encounter_frequencies = iterate_natural_frequencies(FREQ, VEL[0], HEAD[0], A_h, M, C, 9.81, 1e-8)

nat_frequencies = np.power(abs(nat_frequencies_squared), 0.5)


# Compute RAOs
encounter_frequencies, rao = compute_RAOs(VEL, HEAD, FREQ, M, A_h, B_c, B_h, C, F_ex_real, F_ex_imag)
wave_periods = 2*np.pi*np.power(FREQ, -1)
encounter_periods = 2*np.pi*np.power(encounter_frequencies, -1)
# k_encounter = 9.81*np.power(encounter_frequencies, 2)
k = 9.81*np.power(FREQ, 2)


# Plot hydrodynamic coefficients
# Plots added mass in heave
'''
plt.plot(encounter_frequencies, A_h[0, 0, :, 2, 2], 'x-')
plt.xlabel("encounter frequency [rad/s]")
plt.ylabel("$A_{33}^{hyd}$ [kg]")
plt.title("$A_{33}^{hyd}$, " + str(VEL[0]) + " [m/s]")
# plt.vlines([nat_frequencies[0], nat_frequencies[1]], np.amin(A_h[0, 0, :, 2, 2]), np.amax(A_h[0, 0, :, 2, 2]))
plt.grid()
plt.show()
'''
'''
# Plots damping in heave
plt.plot(encounter_frequencies, B_h[0, 0, :, 2, 2], 'x-')
plt.xlabel("encounter frequency [rad/s]")
plt.ylabel("$B_{33}^{hyd}$ [kg/s] ")
plt.title("$B_{33}^{hyd}$, " + str(VEL[0]) + " [m/s]")
# plt.vlines([nat_frequencies[0], nat_frequencies[1]], np.amin(A_h[0, 0, :, 2, 2]), np.amax(A_h[0, 0, :, 2, 2]))
plt.grid()
plt.show()
'''

plt.plot(encounter_periods, A_h[0, 0, :, 4, 4], 'x-')
plt.xlabel("encounter period [s]")
plt.ylabel("$A_{55}^{hyd}$ [kg m^2]")
plt.title("$A_{55}^{hyd}$, " + str(VEL[0]) + " [m/s]")
# plt.vlines([nat_frequencies[0], nat_frequencies[1]], np.amin(A_h[0, 0, :, 2, 2]), np.amax(A_h[0, 0, :, 2, 2]))
plt.grid()
plt.show()

plt.plot(np.divide(encounter_frequencies, 2*np.pi), np.divide(A_h[0, 0, :, 2, 2], 2.4*1025), 'x-')
plt.xlabel("encounter frequency [Hz]")
plt.ylabel("$\\frac{A_{33}^{hyd}}{M^{hyd}}$ [-]")
plt.title("$A_{33}^{hyd}$, " + str(VEL[0]) + " [m/s]")
# plt.vlines([nat_frequencies[0], nat_frequencies[1]], np.amin(A_h[0, 0, :, 2, 2]), np.amax(A_h[0, 0, :, 2, 2]))
plt.grid()
plt.show()


# Plots damping in heave
plt.plot(encounter_frequencies, B_h[0, 0, :, 2, 2], 'x-')
plt.xlabel("encounter frequency [rad/s]")
plt.ylabel("$B_{33}^{hyd}$ [kg/s] ")
plt.title("$B_{33}^{hyd}$, " + str(VEL[0]) + " [m/s]")
# plt.vlines([nat_frequencies[0], nat_frequencies[1]], np.amin(A_h[0, 0, :, 2, 2]), np.amax(A_h[0, 0, :, 2, 2]))
plt.grid()
plt.show()


# Plot RAOs
plot_rao = False
rao_dof = 4  # choose what degree of freedom to plot
DOF_names = ['surge', 'sway', 'heave', 'roll', 'pitch', 'yaw', 'cushion pressure']
xlabel_choice = 3  # chose what to plot the RAO against
xlabel_quantity = ['encounter frequency', 'wave frequency', 'wave periods', 'encounter_periods']

if plot_rao:
    if rao_dof == 0 or rao_dof == 1 or rao_dof == 2:  # Translational degree of freedom
        if xlabel_choice == 0:  # Plots translational DOF against encounter frequency
            plt.plot(encounter_frequencies, abs(rao[rao_dof, :]), 'kx-')
            plt.xlabel("encounter frequency [rad/s]")
        elif xlabel_choice == 1:  # Plots translational DOF against wave frequency
            plt.plot(FREQ, abs(rao[rao_dof, :]), 'kx-')
            plt.xlabel("wave frequency [rad/s]")
        elif xlabel_choice == 2:  # Plots translational DOF against wave periods
            plt.plot(wave_periods, abs(rao[rao_dof, :]), 'kx-')
            plt.xlabel("wave periods [s]")
        elif xlabel_choice == 3:  # Plots translational DOF against encounter periods
            plt.plot(encounter_periods, abs(rao[rao_dof, :]), 'kx-')
            plt.xlabel("encounter periods [s]")
        else:
            raise ValueError

        plt.ylabel('$\eta_{' + str(rao_dof + 1) + '}/ \zeta_a$')  # Sets correct y-label in plot

    elif rao_dof == 3 or rao_dof == 4 or rao_dof == 5:  # Rotational degree of freedom
        if xlabel_choice == 0:  # Plots rotational DOF against encounter frequency
            plt.plot(encounter_frequencies, np.divide(abs(rao[rao_dof, :]), k), 'kx-')
            plt.xlabel("encounter frequency [rad/s]")
        elif xlabel_choice == 1:  # Plots rotational DOF against encounter frequency
            plt.plot(FREQ, np.divide(abs(rao[rao_dof, :]), k), 'kx-')
            plt.xlabel("wave frequency [rad/s]")
        elif xlabel_choice == 2:  # Plots rotational DOF against wave periods
            plt.plot(wave_periods, np.divide(abs(rao[rao_dof, :]), k), 'kx-')
            plt.xlabel("wave periods [s]")
        elif xlabel_choice == 3:  # Plots rotational DOF against encounter periods
            plt.plot(encounter_periods, np.divide(abs(rao[rao_dof, :]), k), 'kx-')
            plt.xlabel("encounter periods [s]")
        else:
            raise ValueError

        plt.ylabel('$\eta_{' + str(rao_dof + 1) + '}/k\zeta_a$')  # Sets correct y-label in plot

    elif rao_dof == 6:  # uniform pressure dof
        if xlabel_choice == 0:  # Plots uniform pressure DOF against encounter frequency
            plt.plot(encounter_frequencies, abs(rao[rao_dof, :]), 'kx-')
            plt.xlabel("encounter frequency [rad/s]")
        elif xlabel_choice == 1:  # Plots uniform pressure DOF against wave frequency
            plt.plot(FREQ, abs(rao[rao_dof, :]), 'kx-')
            plt.xlabel("wave frequency [rad/s]")
        elif xlabel_choice == 2:  # Plots uniform pressure DOF against wave periods
            plt.plot(wave_periods, abs(rao[rao_dof, :]), 'kx-')
            plt.xlabel("wave periods [s]")
        elif xlabel_choice == 3:  # # Plots uniform pressure DOF against encounter periods
            plt.plot(encounter_periods, abs(rao[rao_dof, :]), 'kx-')
            plt.xlabel("encounter periods [s]")
        else:
            raise ValueError
        plt.ylabel('$\eta_{' + str(rao_dof + 1) + '}$')  # Sets correct y-label in plot
    else:
        raise ValueError

    plt.title("RAO in " + DOF_names[rao_dof] + ", " + str(VEL[0]) + " [m/s]")
    plt.grid()
    plt.show()

nat_periods = 2*np.pi*np.power(nat_frequencies, -1)
print(nat_frequencies)
eigen_modes = np.array(eigen_modes)
print(eigen_modes)

print("C_33 = " + str(C_h[0, 0, 50, 2, 2]) + " [N/m]")