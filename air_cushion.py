import numpy as np

def stiffness_matrix_air_cushion(S_0c, h, x_c, z_c, Q_0, dQdp_0, p_0, rho=1025, g=9.81):
    """
    Creates and returns the stiffness matrix containing all terms arising because of the air cushion

    :param S_0c: (double)
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

    # Calculate and add terms to stiffness matrix. (See power point)
    C_c[6, 6] = 0.5*Q_0 - p_0*dQdp_0
    C_c[4, 4] = rho * g * h * S_0c * z_c
    C_c[2, 6] = -rho * g * h * S_0c
    C_c[4, 6] = rho * g * h * S_0c * x_c

    return C_c


def damping_matrix_air_cushion(S_0c, x_c, h, p_0, p_a=101325, gamma=1.4):
    # Initialize damping matrix due to air cushion
    B_c = np.zeros([7, 7])

    # Calculate and add terms to damping matrix. (See power point)
    B_c[6, 6] = p_0*h*S_0c/gamma/(p_0 + p_a)
    B_c[6, 4] = S_0c*x_c  # TODO: Find out what this should be
    B_c[6, 2] = S_0c

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
    S_0c: double
        Total area of the air cushion
    x_c: double
        Distance from centroid to AP

    """

    A_rect = b_c * l_rect  # [m^2] area of rectangle
    A_tri = 0.5 * b_c * l_tri  # [m^2] area of triangle

    # computes area of the air cushion
    S_0c = A_rect + A_tri  # [m^2] total area of air cushion

    # computes distance from AP to centroid
    x_c = (A_rect*0.5*l_rect + A_tri*(l_rect + l_tri/3))/S_0c

    return S_0c, x_c


def interpolate_fan_characteristics(p_0, p, Q):
    # TODO: Implement

    Q_0 = np.interp(p_0, p, Q)  # Interpolates volume flow

    dQdp_0 = 10

    return Q_0, dQdp_0


if __name__ == "__main__":

    # TODO: need to find what these actually are
    l_rect = 9  # [m] length of the rectangular part of the air cushion
    l_tri = 10  # [m] length of the triangular part of the air cushion
    b_c = 4  # [m] beam of the air cushion

    h = 0.4  # [m] mean height between waterline and hull inside air cushion
    z_c = 0.5*h  # [m] vertical centroid of the air cushion relative to the ShipX coordinate system

    p_0 = 3500  # [Pa] excess pressure in air cushion at equilibrium

    # computes air cushion area
    S_0c, x_c = air_cushion_area(l_rect, l_tri, b_c)

    print('Total air cushion area is', S_0c, '[m^2].')
    print('The centroid of the air cushion is located', x_c, '[m] in front of AP.')

    p = np.array([1, 2, 3, 4, 5])
    Q = np.array([1, 2, 3, 4, 5])

    Q_0, dQdp_0 = interpolate_fan_characteristics(p_0, p, Q)

    # computes stiffness matrix
    C_c = stiffness_matrix_air_cushion(S_0c, h, x_c, z_c, Q_0, dQdp_0, p_0)
    print('Stiffness matrix:')
    print(C_c)

    # computes damping matrix
    B_c = damping_matrix_air_cushion(S_0c, x_c, h, p_0)
    print('Damping matrix:')
    print(B_c)
