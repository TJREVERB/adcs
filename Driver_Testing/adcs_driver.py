# Main ADCS Driver
# Python Methods Used:
from getDCM import getDCM
from jd2dvec import jd2dvec
from kep2cart import kep2cart
from kepel import kepel
from q2dcm import q2dcm
from sun_vec import sun_vec
from sunsensors import sunsensors
from utc2jul import utc2jul
from wrldmagm import wrldmagm

import time
import numpy as np
from math import floor, pi, sin, cos, pi, radians, norm
from datetime import datetime, utcnow, date
from orbital import earth, KeplerianElements, utilities
from orbital.utilities import Position, Velocity
from pymap3d import ecef2eci

import logging
from threading import Thread

def main():
    epoch = datetime.utcnow()

if __name__ == "__main__":
    t1 = Thread(target=main, args=(), daemon=True)
    t1.start()