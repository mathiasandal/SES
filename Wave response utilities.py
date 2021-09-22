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


if __name__ == "__main__":

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
    print('magnitude of v4 is ', np.linalg.norm(eigenvectors[:, 3]))