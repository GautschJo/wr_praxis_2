import numpy as np
import linalg

def probe():
    A = np.ones(3,3)
    A = A*5
    print(np.linalg.cholesky(A))