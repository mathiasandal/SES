import numpy as np


def create_mass_matrix(total_mass, r44, r55, r66, r46, x_G, z_G):
    """
    Generates mass matrix for a body with lateral symmetry for all six rigid body motions

    The typical values for radii of gyration in the three rotational degrees of freedom are found from VERES_manual.pdf
    p. 43.

    :param total_mass: double
        total mass of the SES-X vessel without added mass
    :param r44: double
        radii of gyration in roll, typically 0.30*B - 0.45*B
    :param r55: double
        radii of gyration in pitch, typically 0.20*L_pp - 0.30*L_pp
    :param r66: double
        radii of gyration in yaw, typically 0.25*L_pp - 0.30*L_pp
    :param r46: double
        radii of gyration for inertial coupling of yaw and roll, typically zero
    :param: x_G double
        longitudinal position of COG relative to the motion coordinate system
    :param: z_G double
        vertical position of COG relative to the motion coordinate system
    :return: M (6x6) numpy array
        mass matrix containing elements for surge, sway, heave, roll, pitch and yaw
    """

    # Creates mass matrix
    M = total_mass * np.matrix([[1,    0,    0,       0,    z_G,       0],
                                [0,    1,    0,    -z_G,      0,     x_G],
                                [0,    0,    1,       0,   -x_G,       0],
                                [0, -z_G,    0,  r44**2,      0, -r46**2],
                                [z_G,  0, -x_G,       0, r55**2,       0],
                                [0,  x_G,    0, -r46**2,      0,  r66**2]])

    return M


if __name__ == "__main__":
    # location of motion coordinate system relative to intersection of AP, CL and BL
    x_prime = 10  # longitudinal distance from AP to the origin of the motion coordinate system
    z_prime = 3   # vertical distance from BL to the origin of the motion coordinate system

    # main dimensions of BBGreen
    B = 6  # [m] beam of BBGreen
    Lpp = 19.2  # [m] L_pp of BBGreen
    total_mass = 25.6e3  # [kg] total mass of the vessel
    r44 = 0.35*B  # [m] radii of gyration in roll
    r55 = 0.25*Lpp  # [m] radii of gyration in pitch
    r66 = 0.27*Lpp  # [m] radii of gyration in yaw
    r46 = 0  # [m] radii of gyration for inertial coupling of yaw and roll
    lcg = 7.83  # [m] longitudinal center of gravity relative to AP
    vcg = 1.98  # [m] vertical center of gravity relative to BL
    x_G = x_prime - lcg  # [m] longitudinal position of COG relative to the motion coordinate system
    z_G = z_prime - vcg  # [m] vertical position of COG relative to the motion coordinate system

    # Creates mass matrix
    M = create_mass_matrix(total_mass, r44, r55, r66, r46, x_G, z_G)

    print("Mass matrix:")
    print(M)
