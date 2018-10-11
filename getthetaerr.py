# getthetaerr
# Calculates error theta for each axis based on error quaternion
# Translated by Bharath Dileepkumar 10/11/2018
import numpy as np
def getthetaerr(q):
    q = np.transpose(q)
