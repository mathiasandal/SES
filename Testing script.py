from veres import read_re7_file, read_re8_file
from Wave_response_utilities import solve_eq_motion_steady_state, decouple_matrix

import matplotlib.pyplot as plt
import scipy.linalg as la
import numpy as np

# Read *re7 file
VMAS, ADDMAS, DAMP, REST, VEL1, HEAD1, FREQ1 = read_re7_file("Input files/test_input.re7")

# Read *re8 file
REFORCE, IMFORCE, VE2L, HEAD2, FREQ2 = read_re8_file('Input files/test_input.re8')

FORCE = REFORCE + 1j*IMFORCE

#print(REFORCE[:, 0, 0, 0] + 1j*IMFORCE[:, 0, 0, 0])
#print(FORCE[:, 0, 0, 0])

# Gets the complex amplitude of the response eta_j, j = 1,...,6
eta = solve_eq_motion_steady_state(VMAS + ADDMAS[0, 0, 0, :, :], DAMP[0, 0, 0, :, :], REST[0, 0, 0, :, :], FORCE[:, 0, 0, 0], FREQ1[0])

# Frequency index
freq_index = 14

K = REST[0, 0, freq_index, :, :]
A = ADDMAS[0, 0, freq_index, :, :]
M = VMAS + A
B = DAMP[0, 0, freq_index, :, :]
# Calculates eigenvalues and eigenmodes
D, V = la.eig(K, M)  # The warning 'Too many values to unpack do not affect D and V'

print(D)

# Decouple the matrices as surge, heave and pitch is decoupled from sway, roll and yaw due to symmetry
K_surge_heave_pitch = decouple_matrix(K, [0, 2, 4])
M_surge_heave_pitch = decouple_matrix(M, [0, 2, 4])

# Decouple the matrices for heave and pitch
K_heave_pitch = decouple_matrix(K, [2, 4])
M_heave_pitch = decouple_matrix(M, [2, 4])
B_heave_pitch = decouple_matrix(B, [2, 4])

D_surge_heave_pitch, V_surge_heave_pitch = la.eig(K_surge_heave_pitch, M_surge_heave_pitch)

D_heave_pitch, V_heave_pitch = la.eig(K_heave_pitch, M_heave_pitch)

print('For all degrees of freedom:')
print(np.round(D, 1))
print()
print(np.round(V, 1))
print('For surge, heave and pitch only:')
print(np.round(D_surge_heave_pitch, 1))
print()
print(np.round(V_surge_heave_pitch, 1))
print('For heave and pitch only:')
print(np.round(D_heave_pitch, 1))
print()
print(np.round(V_heave_pitch, 1))

print('K:')
print(K_surge_heave_pitch)
print('M:')
print(M_surge_heave_pitch)
print('B:')
print(B_heave_pitch)
print('Natural frequency in heave = ', np.sqrt(K[2, 2]/M[2, 2]))
print('Frequency:', FREQ1[freq_index], 'rad/s')

'''
print(np.sqrt(D[2]))
print()
print(np.round(V, 2))
print()
print(np.round(np.abs(V), 2))
'''

'''
# Plotting coefficients
plt.plot(FREQ1, REST[0, 0, :, 2, 2], 'x-')
plt.xlabel("omega rad/s")
plt.ylabel("C_33")
plt.grid()
plt.show()

plt.plot(FREQ1, DAMP[0, 0, :, 2, 2], 'x-')
plt.xlabel("wave frequency [rad/s]")
plt.ylabel("B_33")
plt.grid()
plt.show()

plt.plot(FREQ1, REFORCE[2, :, 0, 0], 'x-')
plt.xlabel("wave frequency [rad/s]")
plt.ylabel("Force in heave [N/m]")
plt.grid()
plt.show()
'''