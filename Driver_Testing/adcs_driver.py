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
# import gps

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
    return False

def gps_get_cart():
    r = np.array([])
    v = np.array([])
    return np.array([[r],[v]])

def main():
    global epoch
    if(gps_is_on())
        cart = gps_get_cart()
        koe_array = cart2kep(cart[0], cart[1])
    else
        config = load_config('config_adcs.yaml')
        epoch = datetime.utcnow()
        koe_array = np.array([])
        for key, val in config['adcs']['koe'].items():
            koe_array = np.append(koe_array, val)
    print(koe_array)
    config['adcs']['sc']['jd0'] = utc2jul(epoch) # Use write config.
    """
    koe = np.array([])
    for i in range(1,6):
        print(config['adcs']['koe'][i])
        np.append(koe, config['adcs']['koe'][i])
    print(koe)
    """
if __name__ == "__main__":
    t1 = threading.Thread(target=main, args=(), daemon=True)
    t1.start()
    t1.join()
    print("Calculation complete for " + epoch.strftime("%Y-%m-%d %H:%M:%S") + " UTC.")
