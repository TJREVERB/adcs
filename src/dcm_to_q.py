import math
import numpy as np
def dcm2q(dcm):
    q4 = .5*math.sqrt(1+dcm[0][0]+dcm[1][1]+dcm[2][2])
    q1 = (1/(4*q4))*(dcm[1][2]-dcm[2][1])
    q2 = (1/(4*q4))*(dcm[2][0]-dcm[0][2])
    q3 = (1/(4*q4))*(dcm[0][1]-dcm[1][0])
    q = np.array([q1,q2,q3,q4])
    return (q)
