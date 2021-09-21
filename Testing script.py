import numpy as np

A = np.array([[1, 5], [2, 3]])

print(A)
print()
print(np.linalg.norm(np.transpose(A), ord=2, axis=1))
print()
print(np.power(np.linalg.norm(np.transpose(A), ord=2, axis=1), -1))
print()
print(np.diag(np.power(np.linalg.norm(np.transpose(A), ord=2, axis=1), -1)))