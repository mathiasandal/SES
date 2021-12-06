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

    # Reads *.re7 file
    VMAS, ADDMAS, DAMP, REST, VEL_re7, HEAD_re7, FREQ_re7, XMTN_re7, ZMTN_re7 = read_re7_file(path + '//input.re7')

    # Reads *.re8 file
    REFORCE, IMFORCE, VEL_re8, HEAD_re8, FREQ_re8, XMTN_re8, ZMTN_re8 = read_re8_file(path + '//input.re8')

    # Store the necessary data
    A_h = ADDMAS
    B_h = DAMP
    C_h = REST
    F_ex_real = REFORCE
    F_ex_im = IMFORCE

    # Input handling. Checks that the hydrodynamic coefficients from *.re7 have the same velocity, heading and
    # frequencies as the excitation forces in *.re8 file.
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


def iterate_natural_frequencies(wave_frequencies, velocity, heading, added_mass, mass, restoring, g=9.81, tolerance=1e-5):

    """
    Iterates to the correct natural undamped natural frequencies of a system with n degrees of freedom for a SES-X vessel.

    :param wave_frequencies:
    :param velocity:
    :param heading:
    :param added_mass: (NVEL x NHEAD x NFREQ x 6 x 6) Matrix
        Added mass of from the hull calculated using VERES.
    :param mass: (7x7) or (3x3) matrix
        if (7x7)
            contains mass properties of all 7dof, i.e. surge, sway, heave, roll, pitch, yaw and uniform cushion pressure
        if (3x3)
            contains mass properties of three degrees of freedom. Heave, pitch and uniform cushion pressure.
    :param restoring: (3x3) matrix
        if (3x3)
            contains restoring properties of three degrees of freedom. Heave, pitch and uniform cushion pressure.
    :param g: double
        accelerations of gravity
    :param tolerance:

    :return:
    """

    n = len(wave_frequencies)
    m = len(restoring)

    nat_frequencies = np.ones([m], dtype=complex)
    eigen_modes = np.zeros([m, m], dtype=complex)

    # calculates encounter frequency corresponding to each wave frequency and the vessel velocity and wave heading
    encounter_frequencies = wave_frequencies + velocity / g * np.cos(np.deg2rad(heading)) * np.power(wave_frequencies, 2)

    for i in range(m):
        counter = 0  # initialize
        omega_real = 0  # initialize
        index_frequency_upper = int(np.ceil(n/2))
        index_frequency_lower = index_frequency_upper - 1
        err = -1.

        # Get correct size of

        while (counter <= 100 and err >= tolerance) or err == -1:
            # Create matrices for the current encounter frequency
            '''
            M = decouple_matrix(mass, [2, 4, 6])
            A = decouple_matrix(add_row_and_column(added_mass[0, 0, index_frequency, :, :]), [2, 4, 6])
            C = decouple_matrix(restoring, [2, 4, 6])
            '''
            M = mass

            if m == 3:
                if err == -1.:
                    A = decouple_matrix(add_row_and_column(added_mass[0, 0, index_frequency_upper, :, :]), [2, 4, 6])
                else:
                    A_l = decouple_matrix(add_row_and_column(added_mass[0, 0, index_frequency_lower, :, :]), [2, 4, 6])
                    A_u = decouple_matrix(add_row_and_column(added_mass[0, 0, index_frequency_upper, :, :]), [2, 4, 6])
                    A = interpolate_matrices(omega_real, encounter_frequencies[index_frequency_lower], encounter_frequencies[index_frequency_upper], A_l, A_u)

            else:
                if err == -1.:
                    A = add_row_and_column(added_mass[0, 0, index_frequency_upper, :, :])
                else:
                    A_l = add_row_and_column(added_mass[0, 0, index_frequency_lower, :, :])
                    A_u = add_row_and_column(added_mass[0, 0, index_frequency_upper, :, :])
                    A = interpolate_matrices(omega_real, encounter_frequencies[index_frequency_lower], encounter_frequencies[index_frequency_upper], A_l, A_u)

            C = restoring

            nat_freq_temp, eigen_modes_temp = la.eig(C, M + A)

            '''
            if nat_freq_temp[i].real <= 0:
                raise ValueError
            '''

            omega_real = np.sqrt(nat_freq_temp[i].real**2 + nat_freq_temp[i].imag**2)

            # Finds next guess
            dummy_omega, index_frequency = find_closest_value(encounter_frequencies, omega_real)

            if dummy_omega < omega_real:
                index_frequency_lower = index_frequency
                index_frequency_upper = index_frequency + 1
            elif dummy_omega > omega_real:
                index_frequency_lower = index_frequency - 1
                index_frequency_upper = index_frequency

            # Computes relative error
            err = abs((np.sqrt(nat_frequencies[i].real**2 + nat_frequencies[i].imag**2) - omega_real) /
                      np.sqrt(nat_frequencies[i].real**2 + nat_frequencies[i].imag**2))

            # err = abs((encounter_frequencies[index_frequency] - omega_real)/encounter_frequencies[index_frequency])

            # Assigns new natural frequency to array
            nat_frequencies[i] = nat_freq_temp[i]
            eigen_modes[i, :] = eigen_modes_temp[i, :]

            if omega_real == float('inf'):  # infinite frequency
                break  # Breaks if the natural frequency has gone to infinity
            elif index_frequency_lower < 0:
                break  # Breaks if natural frequency goes to zero

            counter += 1  # increment counter

    return nat_frequencies, eigen_modes, encounter_frequencies


def interpolate_matrices(omega, omega_lower, omega_upper, mat_lower, mat_upper):
    """
    Interpolates two (nxn) numpy arrays at two different frequencies.

    :param omega: double
        Frequency of the evaluated frequency
    :param omega_lower: double
        Frequency at the lower end for the interpolation
    :param omega_upper: double
        Frequency at the upper end for the interpolation
    :param mat_lower: (nxn) matrix
        Matrix evaluated at the lower end frequency
    :param mat_upper: (nxn) matrix
        Matrix evaluated at the upper end frequency
    :return:
        mat: (nxn) matrix
            interpolated matrix evaluated at omega.
    """

    # input handling
    if len(mat_upper) != len(mat_lower):
        raise ValueError
    elif omega_lower > omega_upper:
        raise ValueError
    elif omega < omega_lower or omega > omega_upper:
        raise ValueError

    # interpolate
    mat = mat_lower + (omega - omega_lower)/(omega_upper - omega_lower)*(mat_upper - mat_lower)

    return mat


def compute_RAOs(velocity, heading, wave_frequencies, M, A, B_c, B_h, C, F_ex_real, F_ex_imag, f_ex_7, force_combination=1, g=9.81):

    # TODO: Add documentation

    n = len(wave_frequencies)

    encounter_frequencies = wave_frequencies + velocity / g * np.cos(np.deg2rad(heading)) * np.power(wave_frequencies, 2)

    rao = np.zeros([7, n], dtype=complex)

    #force_combination = 2  # 1: Hydro, 2: Wave pumping, 3: both

    for i in range(n):
        A_temp = add_row_and_column(A[0, 0, i, :, :])
        B_temp = add_row_and_column(B_h[0, 0, i, :, :]) + B_c
        # adds a empty column to model no excitation in the air cushion

        if force_combination == 1:
            F_ex_real_temp = np.r_[F_ex_real[:, i, 0, 0], np.array([0])]
            F_ex_imag_temp = np.r_[F_ex_imag[:, i, 0, 0], np.array([0])]
        elif force_combination == 2:
            F_ex_real_temp = np.zeros([7])
            F_ex_imag_temp = np.zeros([7])
            F_ex_real_temp[6] = f_ex_7[i].real
            F_ex_imag_temp[6] = f_ex_7[i].imag
        elif force_combination == 3:
            F_ex_real_temp = np.r_[F_ex_real[:, i, 0, 0], np.array([f_ex_7[i].real])]
            F_ex_imag_temp = np.r_[F_ex_imag[:, i, 0, 0], np.array([f_ex_7[i].imag])]

        # Solves linear system of equations
        rao[:, i] = np.linalg.solve(-encounter_frequencies[i] ** 2 * (M + A_temp) + 1j * encounter_frequencies[i] * B_temp + C, F_ex_real_temp + F_ex_imag_temp*1j)

    return encounter_frequencies, rao


if __name__ == "__main__":
    from Wave_response_utilities import decouple_matrix, add_row_and_column
    from mass_matrix import create_mass_matrix
    from air_cushion import *

    path_veres = 'Input files//Veres input files//22kn//1-15s periods'

    # Read input from .re7 an .re8 files
    A_h, B_h, C_h, F_ex_real, F_ex_im, VEL, HEAD, FREQ, XMTN, ZMTN = read_veres_input(path_veres)

    # Read fan characteristics
    Q, P, rpm = read_fan_characteristics('Input files//fan characteristics//fan characteristics.csv', '1800rpm')

    # Air cushion input variables
    l_rect = 12  # [m] length of the rectangular part of the air cushion
    l_tri = 6  # [m] length of the triangular part of the air cushion
    b_c = 3.4  # [m] beam of the air cushion

    h_b = 0.64  # [m] Cushion plenum height
    h = 0.64  # [m] mean height between waterline(baseline) and hull inside air cushion at AP
    z_c = 0.5 * h  # [m] vertical centroid of the air cushion relative to the ShipX coordinate system

    p_0 = 3500  # [Pa] excess pressure in air cushion at equilibrium

    A_b, x_c = air_cushion_area(l_rect, l_tri, b_c)

    Q_0, dQdp_0 = interpolate_fan_characteristics(p_0, P, Q)



    # main dimensions of BBGreen
    beam = 6  # [m] beam of BBGreen
    Lpp = 19.9  # [m] L_pp of BBGreen

    # location of motion coordinate system relative to intersection of AP, CL and BL
    x_prime = Lpp / 2 + XMTN[0]  # longitudinal distance from AP to the origin of the motion coordinate system
    z_prime = ZMTN[0]  # vertical distance from BL to the origin of the motion coordinate system

    # mass properties
    total_mass = 25.6e3  # [kg] total mass of the vessel
    r44 = 0.35 * beam  # [m] radii of gyration in roll
    r55 = 0.25 * Lpp  # [m] radii of gyration in pitch
    r66 = 0.27 * Lpp  # [m] radii of gyration in yaw
    r46 = 0  # [m] radii of gyration for inertial coupling of yaw and roll
    lcg = 7.83  # [m] longitudinal center of gravity relative to AP
    vcg = 1.98  # [m] vertical center of gravity relative to BL
    x_G = x_prime - lcg  # [m] longitudinal position of COG relative to the motion coordinate system
    z_G = z_prime - vcg  # [m] vertical position of COG relative to the motion coordinate system

    # Create damping and stiffness matrix from air cushion
    B_c = damping_matrix_air_cushion(A_b, x_c, x_prime, h_b, p_0)  # Damping
    C_c = stiffness_matrix_air_cushion(A_b, h, x_c, z_c, x_prime, z_prime, Q_0, dQdp_0, p_0)  # Stiffness

    # Creates mass matrix
    M = create_mass_matrix(total_mass, r44, r55, r66, r46, x_G, z_G)

    # Adds an empty row and column for the uniform pressure degree of freedom eta_7
    M = add_row_and_column(M)

    # Collect stiffness matrices
    C = add_row_and_column(C_h[0, 0, 0, :, :]) + C_c

    decouple = True
    if decouple:
        M = decouple_matrix(M, [2, 4, 6])
        C = decouple_matrix(C, [2, 4, 6])

    nat_frequencies_squared, eigen_modes, encounter_frequencies = iterate_natural_frequencies(FREQ, VEL[0], HEAD[0], A_h, M, C)

    nat_frequencies = np.power(abs(nat_frequencies_squared), 0.5)

    print('Natural frequencies (omega^2):')
    print(nat_frequencies)
    print('Eigen modes:')
    print(eigen_modes)

    plot_added_mass = False

    if plot_added_mass:
        # Plotting coefficients
        plt.plot(FREQ, C_h[0, 0, :, 2, 2], 'x-')
        plt.xlabel("omega rad/s")
        plt.ylabel("C_33")
        plt.title("C_33")
        plt.grid()
        plt.show()

        plt.plot(FREQ, B_h[0, 0, :, 2, 2], 'x-')
        plt.xlabel("wave frequency [rad/s]")
        plt.ylabel("B_33")
        plt.title("B_33, 22kn")
        plt.grid()
        plt.show()

        plt.plot(FREQ, A_h[0, 0, :, 2, 2], 'x-')
        plt.xlabel("wave frequency [rad/s]")
        plt.ylabel("A_33[kg]")
        plt.title("A_33, 22kn")
        plt.grid()
        plt.show()

        plt.plot(encounter_frequencies, A_h[0, 0, :, 2, 2], 'x-')
        plt.xlabel("encounter frequency [rad/s]")
        plt.ylabel("A_33[kg]")
        plt.title("A_33, 22kn")
        plt.vlines([nat_frequencies[0], nat_frequencies[1]], np.amin(A_h[0, 0, :, 2, 2]), np.amax(A_h[0, 0, :, 2, 2]))
        plt.grid()
        plt.show()

        plt.plot(encounter_frequencies, A_h[0, 0, :, 4, 4], 'x-')
        plt.xlabel("encounter frequency [rad/s]")
        plt.ylabel("A_55[kg]")
        plt.title("A_55, 22kn")
        plt.vlines([nat_frequencies[0], nat_frequencies[1]], np.amin(A_h[0, 0, :, 4, 4]),
                   np.amax(A_h[0, 0, :, 4, 4]))
        plt.grid()
        plt.show()

    print('Inertia terms:')
    print('M\t=\t', M[0, 0], '[kg]')
    print('I_55\t=\t', M[1, 1], '[kg m^2]')

    print('Added mass terms:')
    A_33_simplified_1 = np.interp(nat_frequencies[0], encounter_frequencies, A_h[0, 0, :, 2, 2])
    A_33_simplified_2 = np.interp(nat_frequencies[1], encounter_frequencies, A_h[0, 0, :, 2, 2])
    A_55_simplified_1 = np.interp(nat_frequencies[0], encounter_frequencies, A_h[0, 0, :, 4, 4])
    A_55_simplified_2 = np.interp(nat_frequencies[1], encounter_frequencies, A_h[0, 0, :, 4, 4])

    print('A_33 nr1\t=\t', A_33_simplified_1, '[kg]')
    print('A_33 nr2\t=\t', A_33_simplified_2, '[kg]')
    print('A_55 nr1\t=\t', A_55_simplified_1, '[kg m^2]')
    print('A_55 nr2\t=\t', A_55_simplified_2, '[kg m^2]')

    print('Stiffness terms:')

    print('C_33\t=\t', C[0, 0], '[N/m]')
    print('C_55\t=\t', C[1, 1], '[N m^2]')

    print('Uncoupled natural frequencies')
    omega_33_1 = np.sqrt(C[0, 0] / (M[0, 0] + A_33_simplified_1))
    omega_33_2 = np.sqrt(C[0, 0] / (M[0, 0] + A_33_simplified_2))
    omega_55_1 = np.sqrt(C[1, 1] / (M[1, 1] + A_55_simplified_1))
    omega_55_2 = np.sqrt(C[1, 1] / (M[1, 1] + A_55_simplified_2))
    print('Omega_33 nr1\t=\t', omega_33_1, '[rad/s]')
    print('Omega_33 nr2\t=\t', omega_33_2, '[rad/s]')
    print('Omega_55 nr1\t=\t', omega_55_1, '[rad/s]')
    print('Omega_55 nr2\t=\t', omega_55_2, '[rad/s]')

    print(omega_33_1 + omega_33_1**2 * VEL[0]/9.81)
