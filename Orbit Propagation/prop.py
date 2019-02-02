from pyorbital import tlefile
from pyorbital.orbital import Orbital
from pyorbital.orbital import OrbitElements
from datetime import datetime


def runProp():
    #orb = Orbital(tle)
    orb = Orbital("ISS (ZARYA)", tle_file="FILE PATH TO TLE")
    now = datetime.utcnow()
    #print(tle.inclination)
    #print(orb.get_position(now))
    print(orb.get_lonlatalt(now))
    print()
    print(orb.get_next_passes(now, 12, -77.10428, 8.88101, 276, tol=0.001, horizon=0))
if __name__ == '__main__':
    runProp()
