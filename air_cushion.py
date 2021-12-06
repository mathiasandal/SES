import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def stiffness_matrix_air_cushion(A_b, h, x_c, z_c, x_prime, z_prime, Q_0, dQdp_0, p_0, rho=1025, g=9.81):
    """
    Creates and returns the stiffness matrix containing all terms arising because of the air cushion

    :param A_b: (double)
        Total area of the air cushion
    :param h: (double)
        Mean height between waterline and hull inside air cushion
    :param x_c: (double)
        Distance from centroid to AP
    :param z_c: (double)
        Distance from centroid to baseline
    :param Q_0: (double)
        Volume flow of fan at equilibrium
    :param dQdp_0: (double)
        Slope of the fan characteristics curve in the vicinity of the equilibrium volume flow and pressure
    :param p_0: (double)
        Pressure in the air cushion at the equilibrium
    :param rho: (double)
        Density of water
    :param g: (double)
        Acceleration of gravity
    :return:
    C_c: (7x7) numpy array
        Stiffness matrix containing all terms from air cushion
    """
    # Initialize stiffness matrix due to air cushion
    C_c = np.zeros([7, 7])

    # Calculate and add terms to stiffness matrix.
    C_c[6, 6] = 0.5 * Q_0 - p_0 * dQdp_0  # in uniform pressure due to uniform pressure change
    C_c[3, 3] = rho * g * h * (z_prime - z_c) * A_b  # in roll due to roll motion
    C_c[4, 4] = rho * g * h * (z_prime - z_c) * A_b  # in pitch due to pitch motion
    C_c[2, 6] = -rho * g * h * A_b  # in heave due to to uniform pressure change
    C_c[4, 6] = rho * g * h * (x_prime - x_c) * A_b  # in pitch due to uniform pressure change

    return C_c


def damping_matrix_air_cushion(A_b, x_c, x_prime, h_b, p_0, p_a=101325, gamma=1.4):
    """
    Creates and returns the damping matrix containing all terms arising because of the air cushion.

    :param A_b: (float)
        Total area of the air cushion
    :param x_c: (float)
        Distance from centroid to AP
    :param x_prime: (float)
        Distance from AP to motion coordinate system
    :param h_b: (float)
        Mean height between waterline and hull inside air cushion
    :param p_0: (float)
        Pressure in the air cushion at the equilibrium
    :param p_a: (float)
        Atmospheric pressure
    :param gamma: (float)
        Specific heat ratio of air
    :return:
    C_c: (7x7) numpy array
        Damping matrix containing all terms from air cushion
    """

    # Initialize damping matrix due to air cushion
    B_c = np.zeros([7, 7])

    # Calculate and add terms to damping matrix. (See power point)
    B_c[6, 6] = p_0 * h_b * A_b / gamma / (p_0 + p_a)
    B_c[6, 4] = A_b * (x_prime - x_c)
    B_c[6, 2] = A_b

    return B_c


def air_cushion_area(l_rect, l_tri, b_c):
    """
    Computes projected area in zy-plane of a simplified air cushion and its centroid located from AP.
    It is assumed to be a rectangle with a triangle at the end.

            |________l_rect_____| |__l_tri_|
    _        ____________________
    |       |                      *
    |       |       Air cushion       *
    b_c     |             x              >
    |       |                         *
    _       |____________________ *

            |-----x_c----|

    Total area: b_c*l_rect + 0.5 * b_c*l_tri

    :param l_rect: double [m]
        length of rectangular part of the air cushion
    :param l_tri: double [m]
        length of the triangular part of the air cushion
    :param b_c: double [m]
        beam of the air cushion
    :return:
    A_b: double
        Total area of the air cushion
    x_c: double
        Distance from centroid to AP
    """

    A_rect = b_c * l_rect  # [m^2] area of rectangle
    A_tri = 0.5 * b_c * l_tri  # [m^2] area of triangle

    # computes area of the air cushion
    S_0c = A_rect + A_tri  # [m^2] total area of air cushion

    # computes distance from AP to centroid
    x_c = (A_rect * 0.5 * l_rect + A_tri * (l_rect + l_tri / 3)) / S_0c

    return S_0c, x_c


def interpolate_fan_characteristics(p_0, p, Q):
    """
    Calculates the volume flow and the slope of the fan characteristics at a given air cushion pressure and a discrete
    array of pressures and volume flows charcterizing the fan.

    :param p_0: double
        Air cushion pressure at equilibrium
    :param p: double
        Array of pressures that the fan can operate at corresponding to a volume flow
    :param Q: double
        Array of volume flows corresponding to a air cushion pressure
    :return:
    Q_0: double
        Volume flow corresponding to the equilibrium
    dQdp_0: double
        Slope of the fan characteristics in the vicinity of the equilibrium
    """
    # TODO: Fails to interpolate the correct Q. Need to understand why.
    # Q_0 = np.interp(p_0, p, Q)  # Interpolates volume flow # This did not seem to work
    Q_0 = float(interp1d(p, Q)(p_0))  # This works at least for the first test

    dQdp_0 = 0  # initialize slope of fan characteristics

    # finds closest value to p_0 and its index in p
    p_0_closest, p_0_closest_index = find_closest_value(p, p_0)

    # Applies a central finite difference to estimate the slope near the working position (p_0, Q_0)
    if p_0_closest == p_0:
        dQdp_0 = (Q[p_0_closest_index + 1] - Q[p_0_closest_index - 1]) / (
                p[p_0_closest_index + 1] - p[p_0_closest_index - 1])
    elif p_0_closest < p_0:
        alpha = (p_0 - p[p_0_closest_index - 1]) / (p[p_0_closest_index] - p[p_0_closest_index - 1])
        dQdp_0 = (1 - alpha) * (Q[p_0_closest_index] - Q[p_0_closest_index - 2]) / (
                p[p_0_closest_index] - p[p_0_closest_index - 2]) + alpha * (
                         Q[p_0_closest_index + 1] - Q[p_0_closest_index - 1]) / (
                         p[p_0_closest_index + 1] - p[p_0_closest_index - 1])
    elif p_0_closest > p_0:
        alpha = (p_0 - p[p_0_closest_index]) / (p[p_0_closest_index + 1] - p[p_0_closest_index])
        dQdp_0 = (1 - alpha) * (Q[p_0_closest_index + 1] - Q[p_0_closest_index - 1]) / (
                p[p_0_closest_index + 1] - p[p_0_closest_index - 1]) + alpha * (
                         Q[p_0_closest_index + 2] - Q[p_0_closest_index]) / (
                         p[p_0_closest_index + 2] - p[p_0_closest_index])

    return Q_0, dQdp_0


def find_closest_value(arr, val):
    """
    Finds the closest value in a sorted array for a chosen value.

    :param arr: (1xn) numpy array
        Array which function is looping over
    :param val: double
        Value that the function will find the closest value to
    :return:
        closest_value: double
            Closest value the the specified input value
        closest_index: int
            Index of the closest value in the array
    """

    n = len(arr)
    index_closest = 0
    diff_closest = abs(100 * arr[-1] + val)

    for i in range(n):

        if abs(arr[i] - val) < diff_closest:
            index_closest = i
            diff_closest = abs(arr[i] - val)

    closest_value = arr[index_closest]

    return closest_value, index_closest


def read_fan_characteristics(filename, rpm='1800rpm'):
    """
    Reads a csv-file containing the fan characteristics and returns it for a give RPM.

    :param filename: str
        Directory of the csv file containing the fan characteristics.
    :param rpm: str
        String specifying which RPM that should be used.
    :return:
    Q: (nx1) numpy array
        Volume flow values for fan characteristic curve
    P: (nx1) numpy array
        Pressure values for fan characteristic curve
    """

    if rpm not in ['1000rpm', '1400rpm', '1600rpm', '1800rpm', '2108rpm']:
        raise TypeError

    # Reads csv-file and saves it as a pandas DataFrame
    df = pd.read_csv(filename)

    Q = df[['x']].to_numpy()  # gets all Q values from the DataFrame

    P = df[[rpm]].to_numpy()  # gets all P values from the DataFrame

    Q_cut = Q[-1]
    # Cuts the arrays where they intersect P = 0 [Pa]. Found by inspecting fan characteristics manually
    if rpm == '1000rpm':
        Q_cut = 12
    elif rpm == '1400rpm':
        Q_cut = 17
    elif rpm == '1600rpm':
        Q_cut = 19
    elif rpm == '1800rpm':
        Q_cut = 21.6
    elif rpm == '2108rpm':
        Q_cut = 25.5

    # Remove interval that is no well defined
    P = P[Q < Q_cut]
    Q = Q[Q < Q_cut]

    return Q, P, rpm


def wave_pumping_air_cushion(b, l_1, l_2, x_prime, beta, omega, g=9.81):

    # TODO: Add documentation

    k = omega**2/g  # calculate wave number

    # Compute analytical solution to integral of spatial variation over the air cushion
    if np.sin(np.deg2rad(beta)) == 0:  # sin(beta) = 0
        I_1 = b/k/np.cos(np.deg2rad(beta)) * (np.exp(-1j*k*x_prime*np.cos(np.deg2rad(beta))))
        I_2 = 0
    elif np.cos(np.deg2rad(beta)) == 0:  # cos(beta) = 0
        I_1 = 2*np.sin(k*b/2*np.sin(np.deg2rad(beta)))/k/np.sin(np.deg2rad(beta))\
              * (x_prime*np.exp(-1j*k*x_prime*np.cos(np.deg2rad(beta))) -
                (x_prime - l_1)*np.exp(-1j*k*(x_prime - l_1)*np.cos(np.deg2rad(beta))))
        I_2 = 2*l_2/k/np.sin(np.deg2rad(beta))*(1 - np.cos(k*b/2 * np.sin(np.deg2rad(beta))))
    else:  # sin(beta) and cos(beta) is not equal zero
        I_1 = 2 * 1j * np.sin(k * b / 2 * np.sin(np.deg2rad(beta))) / (
                k ** 2 * np.sin(np.deg2rad(beta)) * np.cos(np.deg2rad(beta))) * \
                np.exp(-1j * k * x_prime * np.cos(np.deg2rad(beta)))*(1 - np.exp(1j*k*l_1*np.cos(np.deg2rad(beta))))
        I_2 = 4*l_2*np.exp(-1j*k*(x_prime - l_1)*np.cos(np.deg2rad(beta)))/(k*b**2*np.sin(np.deg2rad(beta))**2 -
              4*l_2**2*k*np.cos(np.deg2rad(beta))**2)*(b/2*np.sin(np.deg2rad(beta)) *
              (np.exp(-1j*k*l_2*np.cos(np.deg2rad(beta))) - np.cos(k*b/2*np.sin(np.deg2rad(beta)))) -
              1j*l_2*np.cos(np.deg2rad(beta))*np.sin(k*b/2*np.sin(np.deg2rad(beta))))

    f_7_hat = 1j*omega*(I_1 + I_2)

    return f_7_hat


if __name__ == "__main__":

    l_1 = 12  # [m] length of the rectangular part of the air cushion
    l_2 = 6  # [m] length of the triangular part of the air cushion
    b = 3.4  # [m] beam of the air cushion

    h_b = 0.64  # [m] Cushion plenum height
    h = 0.5  # [m] mean height between waterline and hull inside air cushion
    z_c = 0.5 * h  # [m] vertical centroid of the air cushion relative to the ShipX coordinate system

    p_0 = 3500  # [Pa] excess pressure in air cushion at equilibrium

    # computes air cushion area
    A_b, x_c = air_cushion_area(l_2, l_2, b)

    # Location of motion coordinate system relative to intersection of AP, CL and BL
    x_prime = 8
    z_prime = 3

    print('Total air cushion area is', A_b, '[m^2].')
    print('The centroid of the air cushion is located', x_c, '[m] in front of AP.')

    # Read fan characteristics from at a specified constant RPM
    Q, P, rpm_dummy = read_fan_characteristics('Input files/fan characteristics/fan characteristics.csv')

    # use interpolation function
    Q_0, dQdp_0 = interpolate_fan_characteristics(p_0, P, Q)

    # plotting results
    plt.plot(P, Q, '-x', label='Q(P)')
    plt.plot(P, Q_0 + dQdp_0 * (P - p_0), 'r', label='Numerical tangent')
    plt.plot(p_0, Q_0, 'k*', label='(p_0, Q_0)')
    plt.xlabel('P [Pa]')
    plt.ylabel('Q [m^3/s]')
    plt.legend(loc='upper right')
    plt.title('Fan characteristics at ' + rpm_dummy[0:-3] + ' RPM')
    plt.show()

    print('Numerical result:')
    print('Q_0 \t=\t', Q_0, '[m^3/s]\ndQdp_0 \t=\t', dQdp_0, '[(m^3s^-1)/(Pa)]')

    # computes stiffness matrix
    C_c = stiffness_matrix_air_cushion(A_b, h, x_c, z_c, x_prime, z_prime, Q_0, dQdp_0, p_0)
    print('Stiffness matrix:')
    print(C_c)

    # computes damping matrix
    B_c = damping_matrix_air_cushion(A_b, x_c, x_prime, h, p_0)
    print('Damping matrix:')
    print(B_c)

    # computes wave pumping
    beta = 90  # [deg] heading of incident waves
    omega = 0.2 * np.pi  # [rad/s] incoming frequency of waves in the beta direction.

    f_7_hat = wave_pumping_air_cushion(b, l_1, l_2, x_prime, beta, omega)
    print(f_7_hat)

    rpms = ['1000rpm', '1400rpm', '1600rpm', '1800rpm', '2108rpm']

    for rpm in rpms:
        Q, P, rpm_plot = read_fan_characteristics("Input files/fan characteristics/fan characteristics.csv", rpm)
        plt.plot(Q, P, label=rpm[:-3] + " " + rpm[-3:].upper())

    plt.title('Fan characteristic curves')
    plt.ylabel('$P$ $[Pa]$')
    plt.xlabel('$Q$ $[m^3/s]$')
    plt.legend(loc='upper right')
    # plt.savefig('Results\\fan_characteristics.pdf')
    plt.show()
