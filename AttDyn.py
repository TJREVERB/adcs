# Attitiude Dynamics
# Function that describes the differential equation of q
# Inputs are:
#   t = time
#   Y = input vector in the form [q, w] where q represents quaternion in
#   vehicle frame and w represents angular velocity in vehicle frame
#   sc = struct of spacecraft properties
#   sim = struct of simulation properties(gain, maximum magnetic moment,
#   transformation matrix from body to magnetorquer frame)
#   KOE = struct of spacecraft keplerian orbital elements
import numpy as np

def att_dyn(t, Y, sc, sim, KOE, jd)
    Y = Y([[q, w]])
    np.reshape(Y, Y.length, order = 'C')
