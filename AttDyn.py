"""
Attitiude Dynamics
Function that describes the differential equation of q
Inputs are:
t = time
Y = input vector in the form [q, w] where q represents quaternion in
vehicle frame and w represents angular velocity in vehicle frame
sc = struct of spacecraft properties
sim = struct of simulation properties(gain, maximum magnetic moment,
transformation matrix from body to magnetorquer frame)
KOE = struct of spacecraft keplerian orbital elements
"""
import numpy as np

def att_dyn(t, Y, sc, sim, KOE, jd):
    Y = np.array([Y]) # Instantiates Y as a 2D array
    Y = Y.reshape(-1, 1) # Transposes Y to a column vector
    #q0 = [[]]
    #for val in Y:
    #    q0.append(val)
    q0 = Y[1:4] # Quaternion rotation from inertial to vehicle frame
    w0 = Y[5:7] # Angular velocity vector of satellite in vehicle frame
    dcm = q2dcm(q0) # DCM from inertial to vehicle frame
