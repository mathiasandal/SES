import numpy as np

def stiffness_matrix_air_cushion():
    # Initialize stiffness matrix due to air cushion
    C_c = np.zeros([7, 7])

    return C_c


def damping_matrix_air_cushion():
    # Initialize damping matrix due to air cushion
    B_c = np.zeros([7, 7])

    return B_c


def air_cushion_area(l_rect, l_tri, b_c):
    """
    Computes projected area in zy-plane of a simplified air cushion and its centroid located from AP.
    It is assumed to be a rectangle with a triangle at the end.

            |________l_rect_____| |__l_tri_|
    _        ____________________
    |       |                      *
    |       |                         *
    b_c     |       Air cushion          >
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


if __name__ == "__main__":

    # TODO: need to find what these actually are
    l_rect = 9  # [m] length of the rectangular part of the air cushion
    l_tri = 10  # [m] length of the triangular part of the air cushion
    b_c = 4  # [m] beam of the air cushion

    # computes air cushion area
    S_0c, x_c = air_cushion_area(l_rect, l_tri, b_c)

    print('Total air cushion area is', S_0c, '[m^2].')
    print('The centroid of the air cushion is located', x_c, '[m] in front of AP.')
