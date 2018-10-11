"""
Attitiude Dynamics
Created by Anonto Zaman, translated to Python by Jason Chen
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
from jdcal import gcal2jd, jd2gcal

def AttDyn(t, Y, sc, sim, KOE, jd):
    Y = np.matrix([Y]) # Instantiates Y as a 2D matrix
    Y = Y.reshape(-1, 1) # Transposes Y to a column vector
    Y = Y.getH() # Returns the complex conjugate transpose, now a row vector
    q0 = Y[1:4] # Quaternion rotation from inertial to vehicle frame
    w0 = Y[5:7] # Angular velocity vector of satellite in vehicle frame
    dcm = q2dcm(q0) # DCM from inertial to vehicle frame

    #Magnetic Field Model
    epochvec = jd2dvec(jd)
    cart = kep2cart(KOE)
    cartloc = cart[1:3]

    getMC
