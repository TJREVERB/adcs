"""
Function BDipole(r, jd, v). Computes the geocentric magnetic field
based on a tilted dipole model. The output is in geocentric inertial
coordinates (ECI). This function includes the effect of a dipole
motion on the Earth.
Inputs:
r       (3,:)   Position vector in the ECI frame (km)
jd      (1,:)   Julian days
v       (3,:)   Velocity (km/s)
Outputs:
b       (3,:)   Magnetic field in the ECI frame (T)
bDot    (3,:)   Derivative of b in the ECI frame (T/s)
Reference:
Wertz, J., ed. "Spacecraft Attitude Determination and
Control," Kluwer, 1976, 783-784.
Created by Jason Chen 10/18/19
"""
import numpy as np
import math

def BDipole(*args):
    if(len(args) < 2):
        a = 7000;
        el = [a, math.radians(55), 0, 0, 0, 0, 0]
        p = Period(a)
    #bDot = np.zeros(3,n)
