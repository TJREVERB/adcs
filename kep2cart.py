"""
Function kep2cart(KOE), returns Cartesian position/velocity vectors
[R;V] from Keplerian orbital elements (KOE), using the EGM-96 value
of Earth's gravitational constant, GM = 3.986004415E+14 m^3/sec^2.
KOE is defined as in kepel.py, as follows:
    KOE.sma = semi-major axis
    KOE.ecc = eccentricity
    KOE.incl = inclination
    KOE.raan = right ascension of the ascending node
    KOE.argp = argument of periapse
    KOE.tran = true anomaly
---MORE DESCRIPTION---
Created by Jason Chen 10/6/18
"""
import poliastro
from astropy import coordinates
