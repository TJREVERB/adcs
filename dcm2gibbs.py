# dcm2gibbs
# created by Anonto Zaman, translated to python by Kirthi Kumar
# Conversion of direction cosine matrix to 3-1-2 gibbs vector
# Output is a row vector with components representing rotations about the
# 3-1-2 axes respectively

import numpy as npy
import math

def dcm2gibbs(dcm):
    g = ((1, 3))
    g.zeros(3)
    g[1, 1] = -1 *math.atan(dcm[2,1]/dcm[2,2])
    g[1, 2] = math.asin(dcm[2, 3])
    g[1, 3] = math.atan(dcm[1,3]/dcm[3,3])
    return g