import numpy as np
from scipy.linalg import *

A = np.array([[1,3],[2,2]])

eig_val, eig_vec = np.linalg.eig(A)

print(eig_val)
print(eig_vec[:,0])
