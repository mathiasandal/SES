import numpy as np
import scipy.linalg as la


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
    eta = np.linalg.solve(-omega ** 2 * M + 1j * omega * B + C, F)

    return eta


def calculate_eigenfrequency_and_eigenmodes(M, C):

    """
    Calculates the eigenvalues and the corresponding eigenvectors og the mass and stiffness matrix.

    It is used a Cholesky decomposition, meaning that each eigenvector is orthogonal. To go back to the physical space, x,
    the eigenvector needs to be transformed like this: x = inv(transpose(L))*q, where q is the eigenvector.

    :param M: Mass matrix, nxn numpy.array
    :param C: Stiffness matrix, nxn numpy.array
    :return:
        eigenvalues: 1xn array with eigenvalues
        eigenvectors: nxn array with eigenvectors. Each column vector is the normalized eigenvector corresponding to the
        eigenvalue at the same index in eigenvalues.
    """

    # Cholesky decomposition
    L = np.linalg.cholesky(M)  # L decomposition of M, M = L*transpose(L)
    L_T = np.transpose(L)  # Transpose of L,
    L_inv = np.linalg.inv(L)  # Inverse of L
    L_T_inv = np.linalg.inv(L_T)  # Inverse of transpose(L)

    A = L_inv @ C @ L_T_inv  # Matrix multiplication

    # A = np.matmul(L_inv, np.matmul(C, L_T_inv))

    eigenvalues, eigenvectors_orthogonal = np.linalg.eig(A)  # Solve for eigenvalues of the orthogonal eigenmodes

    # Transform from Cholesky to physical space (x = L_T_inv*q)
    # eigenvectors = np.matmul(L_T_inv, eigenvectors_orthogonal)
    eigenvectors = L_T_inv @ eigenvectors_orthogonal

    # Mulitplying eigenvectors by 1/norm(vec[:, i]) to normalize the eigenvectors
    RHS = np.diag(np.power(np.linalg.norm(np.transpose(eigenvectors), ord=2, axis=1), -1))
    print(RHS)
    eigenvectors_normalized = eigenvectors @ RHS

    ''''''
    # Testing
    eigenvector1 = eigenvectors_normalized[:, 0]
    eigenvector2 = eigenvectors_normalized[:, 1]


    print('Eigenvectors: \n', eigenvectors_normalized)
    print()
    print('Eigenvector1 = ', eigenvector1)
    print('Magnitude of Eigenvector1 = ', np.linalg.norm(eigenvector1, ord=2))
    print('Eigenvector2 = ', eigenvector2)
    print('Magnitude of Eigenvector2 = ', np.linalg.norm(eigenvector2, ord=2))
    print()

    return eigenvalues, eigenvectors_normalized


def decouple_matrix(mat_in, elem):
    """
    Reduce the size of the matrix with the elements specified in elem

    :param mat_in: an (MXM) numpy array
    :param elem: (1XN) numpy array containing the indexes to be put in mat_out
    :return:
        mat_out: (NXN) numpy array containing elements from mat_in specified by elem
    """

    n = len(elem)
    mat_out = np.zeros([n, n])

    for i in range(n):
        for j in range(n):
            mat_out[i, j] = mat_in[elem[i], elem[j]]

    return mat_out


def add_row_and_column(mat):
    """
    Adds a row and column with zeros to a matrix
    :param mat:
    :return:
    """

    n = len(mat)

    return np.r_[np.c_[mat, np.zeros(n)], [np.zeros(n+1)]]


def cushion_terms(gamma, rho, g, h, h_b, p_0, p_a, Q_0, dQdp0, S_0):
    """
    Calculates and returns terms in the coupled equations due to the air cushion.

    :param gamma: [-] ratio of specific heat for air
    :param rho: [kg/m^3] density of water
    :param g: [m/s^2] gravity
    :param h: [m] difference of mean water level inside and outside the air cushion
    :param h_b: [m] distance from water level inside air cushion to the hull
    :param p_0: [Pa] excess pressure inside the air cushion
    :param p_a: [Pa] atmospheric pressure
    :param Q_0: [m^3/s] mean volume flow of the fan
    :param dQdp0: # [(m^3s^-1)/(Pa)] Slope of fan characteristic curve estimated from Figure 5.6 in Faltinsen High-Speed
    :param S_0: # [m^2] Area of the waterline inside the air cushion
    :return:
    """
    # TODO: Should the function take an array of Q and P(fan characteristics) or just get the interpolated value?
    C_cushion = np.zeros([7, 7])  # initialize stiffness matrix to contain cushion terms
    B_cushion = np.zeros([7, 7])  # initialize damping matrix to contain cushion terms

    K_1 = p_0*h_b*S_0/(gamma*(p_0 + p_a))
    K_2 = 0.5*Q_0 - p_0*dQdp0
    A_0 = S_0  # [m^2] mean free surface area
    A_1 = 0  # [m^3](?) longitudinal moment of area of the cushion

    # insert terms in matrices in the continuity of air equation, see p. 21 in VERES SES Manual
    B_cushion[6, 2] = A_0
    B_cushion[6, 4] = -A_1
    B_cushion[6, 6] = K_1
    C_cushion[6, 6] = K_2

    # hydrostatic terms coming from the air cushion, see p. 15 in VERES SES Manual
    # TODO: Investigate these terms
    C_cushion[2, 6] = 0
    C_cushion[4, 6] = 0

    return C_cushion, B_cushion

''''''
def iterate_natural_frequencies(wave_frequencies, velocity, heading, added_mass, mass, restoring, g=9.81):
    n = len(wave_frequencies)
    m = len(mass)
    nat_frequencies = np.zeros([n, m], dtype=complex)
    eigen_modes = np.zeros([n, m, m], dtype=complex)

    # calculates encounter frequency corresponding to each wave frequency and the vessel velocity and wave heading
    encounter_frequencies = wave_frequencies + velocity/g*np.cos(np.deg2rad(heading))*np.power(wave_frequencies, 2)

    index_nat_freq = np.zeros([1, m])  # containing index of natural frequency that is closest to the encounter frequency it is calculated from
    diff_nat_encounter_freq = np.zeros([n, m])

    for i in range(n):
        nat_frequencies[i], eigen_modes[i] = la.eig(restoring[i], mass + added_mass[i])

    return nat_frequencies, eigen_modes, wave_frequencies, encounter_frequencies


#def iterate_natural_frequencies(wave_frequencies, velocity, heading, added_mass, mass, restoring, g=9.81):



if __name__ == "__main__":
    from veres import read_re7_file

    '''
    # Mass matrix
    M = np.array([[1.0000, 0.0000, 0.0000, 0.0000],
                  [0.0000, 2.0000, 0.0000, 0.0000],
                  [0.0000, 0.0000, 8.0000, 0.0000],
                  [0.0000, 0.0000, 0.0000, 4.0000]])

    # Stiffness matrix
    K = np.array([[1.0000, 0.5000, 0.3333, 0.2500],
                  [0.5000, 1.0000, 0.6667, 0.5000],
                  [0.3333, 0.6667, 1.0000, 0.7500],
                  [0.2500, 0.5000, 0.7500, 1.0000]])

    D, V = la.eig(K, M)  # la.eig(K, M) does the exact same as calculate_eigenfrequency_and_eigenmodes(M, C): hahaha

    # Testing calculate_eigenfrequency_and_eigenmodes(M, C)

    eigenvalues, eigenvectors = calculate_eigenfrequency_and_eigenmodes(M, K)

    print('Eigenvalues:')
    print(eigenvalues)
    print()
    print('Eigenvectors:')
    print(eigenvectors)
    print()
    print('D: ')
    print(D)
    print()
    print('V: ')
    print(V)
    print()

    print('v1 = ', eigenvectors[:, 0])
    print('v2 = ', eigenvectors[:, 1])
    print('v3 = ', eigenvectors[:, 2])
    print('v4 = ', eigenvectors[:, 3])

    print('v1 * v2 = ', np.matmul(eigenvectors[:, 0], eigenvectors[:, 1]))
    print('v1 * v3 = ', np.matmul(eigenvectors[:, 0], eigenvectors[:, 2]))
    print('v1 * v4 = ', np.matmul(eigenvectors[:, 0], eigenvectors[:, 3]))
    print('v2 * v3 = ', np.matmul(eigenvectors[:, 1], eigenvectors[:, 2]))
    print('v2 * v4 = ', np.matmul(eigenvectors[:, 1], eigenvectors[:, 3]))
    print('v3 * v4 = ', np.matmul(eigenvectors[:, 2], eigenvectors[:, 3]))
    print()
    print('magnitude of v1 is ', np.linalg.norm(eigenvectors[:, 0]))
    print('magnitude of v2 is ', np.linalg.norm(eigenvectors[:, 1]))
    print('magnitude of v3 is ', np.linalg.norm(eigenvectors[:, 2]))
    print('magnitude of v4 is ', np.linalg.norm(eigenvectors[:, 3]))'''

    # Test iterate_natural_frequencies(omegas, velocity, heading, added_mass, restoring):
    VMAS, ADDMAS, DAMP, REST, VEL, HEAD, FREQ, XMTN, ZMTN = read_re7_file("Input files/test_input.re7")

    # indices
    iheading = 0
    ivel = 0

    nat_frequencies, eigen_modes, wave_frequencies, encounter_frequencies = iterate_natural_frequencies(FREQ, VEL[ivel], HEAD[iheading], ADDMAS[ivel, iheading, :, :, :], VMAS, REST[ivel, iheading, :, :, :])

    print(nat_frequencies[14])
