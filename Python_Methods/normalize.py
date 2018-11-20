import numpy as np
from numpy import linalg

def normalize(array):
    return np.linalg.norm(array, axis=-1)[:, np.newaxis]
"""
array = np.array([[ 0.        ,  0.4472136 ,  0.89442719],
        [ 0.42426407,  0.56568542,  0.70710678],
        [ 0.49153915,  0.57346234,  0.65538554]])
print(array/normalize(array))
"""
