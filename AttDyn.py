"""
Attitiude Dynamics
Created by Anonto Zaman, BEING TRANSLATED by Jason Chen
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
import datetime
from jdcal import gcal2jd, jd2gcal

def att_dyn(t, Y, sc, sim, KOE, jd):
    Y = np.array([Y]) # Instantiates Y as a 2D array
    Y = Y.reshape(-1, 1) # Transposes Y to a column vector
    #q0 = [[]]
    #for val in Y:
    #    q0.append(val)
    q0 = Y[1:4] # Quaternion rotation from inertial to vehicle frame
    w0 = Y[5:7] # Angular velocity vector of satellite in vehicle frame
    dcm = q2dcm(q0) # DCM from inertial to vehicle frame

    #Magnetic Field Model
    ps = jd - 2400000.5 # Done to increase precision to microseconds
    dtarray = [] # Blank datetime array
    datetime.datetime(jd2gcal(2400000.5, ps)))
    dtarray.append
    # No idea if this works, adds converted Julian date to the Gregorian
    # format of (YYYY, M, D, decimal of a day). Instantiates this as a
    # datetime object. Appends this object to a datetime array.
    for x in range(0, len(dtarray))
        epochvec = [[]]
