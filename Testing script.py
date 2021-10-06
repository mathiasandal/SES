import numpy as np
from Wave_response_utilities import decouple_matrix

A = np.array([[1., 5., 8.], [2., 3., 16.], [32., 54., 90.]])
'''
print(A)
print()
print(np.linalg.norm(np.transpose(A), ord=2, axis=1))
print()
print(np.power(np.linalg.norm(np.transpose(A), ord=2, axis=1), -1))
print()
print(np.diag(np.power(np.linalg.norm(np.transpose(A), ord=2, axis=1), -1)))
'''

M = decouple_matrix(A, [1, 2])

print(A)
print()
print(M)