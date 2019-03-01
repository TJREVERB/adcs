from datetime import datetime
import numpy as np
from pyorbital import tlefile
from pyorbital.orbital import Orbital
from pyorbital.orbital import OrbitElements
from datetime import datetime

def get_lla(t):
    orb = Orbital("TJREVERB", tle_file="tjreverbtle.txt")
    return({'lat':orb.get_lonlatalt(t)[0],'lon':orb.get_lonlatalt(t)[1],'alt':orb.get_lonlatalt(t)[2]})
