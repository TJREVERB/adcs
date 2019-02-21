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
from cart2kep import cart2kep
from gps import get_data

import time
import numpy as np
from numpy import linalg
from math import floor, pi, sin, cos, pi, radians
from datetime import datetime, date
from orbital import earth, KeplerianElements, utilities
from orbital.utilities import Position, Velocity
from pymap3d import ecef2eci

import logging
import threading
import yaml

def load_config(config_file):
    with open(config_file, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as error:
            print(error)

def write_config(config_file, data):
    with open(config_file, 'w') as stream:
        yaml.dump(data, stream, default_flow_style=False)

def gps_is_on():
    return True

def tle_get_data():
    return None;

def main():
    global epoch
    # If GPS is on, get Cartesian (position, velocity) vectors and UTC time.
    # Convert Cartesian coordinates and time to a Keplerian Elements array.
    # Update the config_adcs.yaml file with new KOE array.
    if(gps_is_on()):
        data = gps.get_data()
        cart = data['']
        koe_array = cart2kep(cart[0], cart[1])
        epoch = data['time']
    else:
        # write_config('config_adcs.yaml', tle_get_data())
        config = load_config('config_adcs.yaml')
        epoch = datetime.utcnow()
        koe_array = np.array([])
        for key, val in config['adcs']['koe'].items():
            koe_array = np.append(koe_array, val)
    print(koe_array)
    config['adcs']['sc']['jd0'] = utc2jul(epoch) # Use write config.

if __name__ == "__main__":
    t1 = threading.Thread(target=main, args=(), daemon=True)
    t1.start()
    t1.join()
    print("Calculation complete for " + epoch.strftime("%Y-%m-%d %H:%M:%S") + " UTC.")
