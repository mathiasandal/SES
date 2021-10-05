from read_re8_file import read_re8_file
from read_re7_file import read_re7_file
from Wave_response_utilities import solve_eq_motion_steady_state
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

D, V = la.eig(REST[0, 0, 0, :, :], VMAS + ADDMAS[0, 0, 0, :, :])

print(np.round(np.abs(V), 3))
