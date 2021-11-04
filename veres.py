import numpy as np
import scipy.linalg as la
from air_cushion import find_closest_value
from Wave_response_utilities import decouple_matrix, add_row_and_column

'''Contains all functions related to things retrieved from VERES'''


def read_re7_file(filename):
    # TODO: Fix/add documentation
    """
    Reads a *.re7 file created by VERES

    A description of a *.re7 file is found in page 114 of the ShipX Vessel Responses (VERES) User's manual.
    """

    f = open(filename, "r")  # open file for reading

    for i in range(6):  # skips first six lines
        f.readline()

    # Read in parameters
    RHOSW, GRAV = [float(i) for i in f.readline().split()]
    LPP, BREADTH, DRAUGHT = [float(i) for i in f.readline().split()]
    LCG, VCG = [float(i) for i in f.readline().split()]

    # saves line looking like "-1 FILEVER (=2)" in Veres_Manual.pdf
    FILEVER = [int(i) for i in f.readline().split()]

    # Read in number of velocities, headings, frequencies and degrees of freedom
    NOVEL, NOHEAD, NOFREQ, NDOF = [int(i) for i in f.readline().split()]

    # Initialize vectors to contain velocities, headings and frequencies
    VEL = np.zeros(NOVEL)  # [m/s] Velocity of the vessel
    HEAD = np.zeros(NOHEAD)  # [deg] Wave heading
    FREQ = np.zeros(NOFREQ)  # [rad/s] Wave frequency

    # Initialize values to be stores for each velocity
    SINK = np.zeros(NOVEL)  # [m] Sink
    TRIM = np.zeros(NOVEL)  # [deg] Trim
    XMTN = np.zeros(NOVEL)  # [m] X-pos. of the motion coordinate system (relative to Lpp/2)
    ZMTN = np.zeros(NOVEL)  # [m] Z-pos. of the mo􀆟on coordinate system (relative to BL)

    # Read in mass matrix
    VMAS = np.zeros([NDOF, NDOF])  # Initialize mass matrix

    for j in range(NDOF):  # Read and converts each line from string to floats
        VMAS[j, :] = [float(i) for i in f.readline().split()]

    # Have one (NDOF x NDOF) matrix for each combinations of velocity, heading and wave frequency
    ADDMAS = np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Hydrodynamic added mass
    ADDADDMAS = np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Additional added mass
    DAMP = np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Hydrodynamic damping
    ADDDAMP = np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Additional damping
    REST = np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Restoring
    ADDREST = np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Additional restoring

    # Initialize viscous roll damping
    VISCDL = np.zeros([NOVEL, NOHEAD, NOFREQ])  # Viscous roll damping, linear part
    VISCDN = np.zeros([NOVEL, NOHEAD, NOFREQ])  # Viscous roll damping, nonlinear part
    VISCDNL = np.zeros([NOVEL, NOHEAD, NOFREQ])  # Viscous roll damping, linearized

    for i in range(NOVEL):
        VEL[i], SINK[i], TRIM[i], XMTN[i], ZMTN[i] = [float(m) for m in f.readline().split()]
        for j in range(NOHEAD):
            HEAD[j] = float(f.readline())  # Should only contain one element
            for k in range(NOFREQ):
                FREQ[k] = float(f.readline())  # Should only contain one element

                for m in range(NDOF):
                    ADDMAS[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                for m in range(NDOF):
                    ADDADDMAS[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                for m in range(NDOF):
                    DAMP[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                for m in range(NDOF):
                    ADDDAMP[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                for m in range(NDOF):
                    REST[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                for m in range(NDOF):
                    ADDREST[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                # Read in viscous roll damping
                VISCDL[i, j, k], VISCDN[i, j, k], VISCDNL[i, j, k] = [float(i) for i in f.readline().split()]
    return VMAS, ADDMAS, DAMP, REST, VEL, HEAD, FREQ, XMTN, ZMTN


def read_re8_file(filename):
    # TODO: Fix/add documentation
    """
    Reads a *.re8 file created by VERES

    A description of a *.re8 file is found in page 116 of the ShipX Vessel Responses (VERES) User's manual.

    Only works for High-Speed formulation in VERES
    """

    f = open(filename, 'r')

    for i in range(6):  # skips first six lines
        f.readline()

    # Read in parameters
    RHOSW, GRAV = [float(i) for i in f.readline().split()]
    LPP, BREADTH, DRAUGHT = [float(i) for i in f.readline().split()]
    LCG, VCG = [float(i) for i in f.readline().split()]

    # Read in number of velocities, headings, frequencies and degrees of freedom
    NOVEL, NOHEAD, NOFREQ, NDOF = [int(i) for i in f.readline().split()]


    # Initialize vectors to contain velocities, headings and frequencies
    VEL = np.zeros(NOVEL)  # [m/s] Velocity of the vessel
    HEAD = np.zeros(NOHEAD)  # [deg] Wave heading
    FREQ = np.zeros(NOFREQ)  # [rad/s] Wave frequency

    # Initialize values to be stores for each velocity
    SINK = np.zeros(NOVEL)  # [m] Sink
    TRIM = np.zeros(NOVEL)  # [deg] Trim
    XMTN = np.zeros(NOVEL)  # [m] X-pos. of the motion coordinate system (relative to Lpp/2)
    ZMTN = np.zeros(NOVEL)  # [m] Z-pos. of the mo􀆟on coordinate system (relative to BL)

    # Initialize Force components
    # Real parts
    REFORCE = np.zeros([NDOF, NOFREQ, NOHEAD, NOVEL])  # Real part
    IMFORCE = np.zeros([NDOF, NOFREQ, NOHEAD, NOVEL])  # Imaginary part

    for i in range(NOVEL):
        VEL[i], SINK[i], TRIM[i], XMTN[i], ZMTN[i] = [float(m) for m in f.readline().split()]
        for j in range(NOHEAD):
            HEAD[j] = float(f.readline())  # Should only contain one element
            for k in range(NOFREQ):
                FREQ[k] = float(f.readline())  # Should only contain one element

                for m in range(NDOF):
                    REFORCE[m, k, j, i], IMFORCE[m, k, j, i] = [float(m) for m in f.readline().split()][1:]

    return REFORCE, IMFORCE, VEL, HEAD, FREQ, XMTN, ZMTN


def read_veres_input(path):
    # TODO: Add documentation

    VMAS, ADDMAS, DAMP, REST, VEL_re7, HEAD_re7, FREQ_re7, XMTN_re7, ZMTN_re7 = read_re7_file(path + '//input.re7')

    REFORCE, IMFORCE, VEL_re8, HEAD_re8, FREQ_re8, XMTN_re8, ZMTN_re8 = read_re8_file(path + '//input.re8')

    A_h = ADDMAS
    B_h = DAMP
    C_h = REST
    F_ex_real = REFORCE
    F_ex_im = IMFORCE

    if VEL_re7.all() == VEL_re8.all() and HEAD_re7.all() == HEAD_re8.all() and FREQ_re7.all() == FREQ_re8.all() and \
            XMTN_re7.all() == XMTN_re8.all() and ZMTN_re7.all() == ZMTN_re8.all():
        VEL = VEL_re7
        HEAD = HEAD_re7
        FREQ = FREQ_re7
        XMTN = XMTN_re7
        ZMTN = ZMTN_re7
    else:
        raise ValueError

    return A_h, B_h, C_h, F_ex_real, F_ex_im, VEL, HEAD, FREQ, XMTN, ZMTN


def iterate_natural_frequencies(wave_frequencies, velocity, heading, added_mass, mass, restoring, g=9.81, tolerance=1e-6):

    n = len(wave_frequencies)
    m = len(restoring)

    nat_frequencies = np.zeros([m], dtype=complex)
    eigen_modes = np.zeros([m, m], dtype=complex)

    # calculates encounter frequency corresponding to each wave frequency and the vessel velocity and wave heading
    encounter_frequencies = wave_frequencies + velocity / g * np.cos(np.deg2rad(heading)) * np.power(wave_frequencies, 2)

    for i in range(m):
        counter = 0
        index_frequency = int(np.ceil(n/2))
        err = -1.

        # Get correct size of

        while (counter <= 100 and err <= tolerance) or err == -1:
            # Create matrices for the current encounter frequency
            '''
            M = decouple_matrix(mass, [2, 4, 6])
            A = decouple_matrix(add_row_and_column(added_mass[0, 0, index_frequency, :, :]), [2, 4, 6])
            C = decouple_matrix(restoring, [2, 4, 6])
            '''
            M = mass
            A = add_row_and_column(added_mass[0, 0, index_frequency, :, :])
            C = restoring

            nat_freq_temp, eigen_modes_temp = la.eig(C, M + A)

            if nat_freq_temp[i].real <= 0:
                raise ValueError

            # Computes relative error
            err = abs((encounter_frequencies[index_frequency] - np.sqrt(nat_freq_temp[i].real))/encounter_frequencies[index_frequency])

            # Finds next guess
            dummy_frequency, index_frequency = find_closest_value(encounter_frequencies, np.sqrt(nat_freq_temp[i].real))

            # Assigns new natural frequency to array
            nat_frequencies[i] = nat_freq_temp[i]
            eigen_modes[i, :] = eigen_modes_temp[i, :]

            counter += 1  # increment counter

    return nat_frequencies, eigen_modes



if __name__ == "__main__":
    from Wave_response_utilities import decouple_matrix, add_row_and_column
    from mass_matrix import create_mass_matrix
    from air_cushion import *

    path_veres = 'Input files//Veres input files'

    # Read input from .re7 an .re8 files
    A_h, B_h, C_h, F_ex_real, F_ex_im, VEL, HEAD, FREQ, XMTN, ZMTN = read_veres_input(path_veres)

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

    # Create damping and stiffness matrix from air cushion #TODO: input should be difference between VERES coordinate system and cushion centroid
    B_c = damping_matrix_air_cushion(S_0c, x_c - XMTN, h, p_0)  # Damping
    C_c = stiffness_matrix_air_cushion(S_0c, h, x_c - XMTN, z_c - ZMTN, Q_0, dQdp_0, p_0)  # Stiffness

    # main dimensions of BBGreen
    B = 6  # [m] beam of BBGreen
    Lpp = 19.2  # [m] L_pp of BBGreen
    total_mass = 25.6e3  # [kg] total mass of the vessel
    r44 = 0.35*B  # [m] radii of gyration in roll
    r55 = 0.25*Lpp  # [m] radii of gyration in pitch
    r66 = 0.27*Lpp  # [m] radii of gyration in yaw

    # Creates mass matrix
    M = create_mass_matrix(total_mass, r44, r55, r66)

    nat_frequencies, eigen_modes = iterate_natural_frequencies(FREQ, VEL[0], HEAD[0], A_h, M, add_row_and_column(C_h[0, 0, 0, :, :]) + C_c, g=9.81, tolerance=1e-3)

    print(nat_frequencies)
    print(eigen_modes)
    '''
    REFORCE, IMFORCE, VEL, HEAD, FREQ = read_re8_file('Input files/test_input.re8')

    VMAS, ADDMAS, DAMP, REST, VEL7, HEAD7, FREQ7 = read_re7_file("Input files/test_input.re7")

    print(REFORCE[0, 0, 0, 0])
    '''