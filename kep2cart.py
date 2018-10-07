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
"""
TO USE ASTROPY:
pip install -U pytest
pip install clang
pip install --upgrade setuptools
pip install astropy --no-deps
IF YOU ENCOUNTER THIS ERROR: Microsoft Visual C++ 14.0 is required.
GO TO: https://visualstudio.microsoft.com/visual-cpp-build-tools/
INSTALL: Build Tools for Visual Studio 2017
ONCE THE EXE IS RUNNING, INSTALL THE CORE C++ PACKAGE (around 5.3 GB).
TRY AGAIN: pip install astropy
pip install poliastro
"""
import poliastro
import astropy

#Test input: (7712.1869, 0.001, 63.43, 135.00, 90.00, 0.0)
#Test output: (-2436.45, -2436.45, 6891.037)
