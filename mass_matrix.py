import numpy as np


def create_mass_matrix(total_mass, r44, r55, r66):

    '''
    Generates mass matrix for a SES-X vessel including all rigid body motions and uniform pressure

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
    :return: M (7x7) numpy array
        mass matrix containing elements for surge, sway, heave, roll, pitch, yaw and uniform cushion pressure
    '''

    M = np.zeros([7, 7])

    M[0, 0] = total_mass             # M_11, surge
    M[1, 1] = total_mass             # M_22, sway
    M[2, 2] = total_mass             # M_33, heave
    M[3, 3] = total_mass * r44 ** 2  # I_44, roll
    M[4, 4] = total_mass * r55 ** 2  # I_55, pitch
    M[5, 5] = total_mass * r66 ** 2  # I_66, heave

    return M


if __name__ == "__main__":

    B = 6  # [m] beam of BBGreen
    Lpp = 19.2  # [m] L_pp of BBGreen
    total_mass = 25.6e3  # [kg] total mass of the vessel
    r44 = 0.35*B  # [m] radii of gyration in roll
    r55 = 0.25*Lpp  # [m] radii of gyration in pitch
    r66 = 0.27*Lpp  #

    # Creates mass matrix
    M = create_mass_matrix(total_mass, r44, r55, r66)

    print("Mass matrix:")
    print(M)