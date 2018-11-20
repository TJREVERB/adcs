# getthetaerr
# Calculates error theta for each axis based on error quaternion
# Translated by Bharath Dileepkumar 10/18/2018
import numpy as np
def getthetaerr(q):
    q = np.matrix(q) # Initializes q into a row vector (2d matrix in Python)
    thetaerr = 2 * (q[0:3]/q[3]) # Calculates attitude error in each axis by dividing q by scalar component
    thetaerr = thetaerr.transpose()
    thetaerr.shape = (1,3)
    return thetaerr
    
