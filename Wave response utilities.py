import numpy as np


def solve_eq_motion_steady_state(M, B, C, F, omega):
    """
    Solves the steady state solution for a given frequency.

    :param M: Mass matrix, nxn numpy array
    :param B: Damping matrix, nxn numpy array
    :param C: Restoring matrix, nx1 numpy array
    :param F: Excitation amplitudes, nx1 numpy array
    :param omega: [rad/s] excitation frequency
    :return: eta: response amplitudes
    """
    eta = np.linalg.solve(-omega**2 * M + 1j*omega * B + C, F)

    return eta


def calculate_eigenfrequency_and_eigenmodes(M, C):

    omegas = 0

    return omegas