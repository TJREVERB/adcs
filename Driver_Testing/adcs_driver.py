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
from gps_dummy import get_data

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
        try:
            yaml.dump(data, stream, default_flow_style=False)
        except yaml.YAMLError as error:
            print(error)

def gps_is_on():
    return True

def tle_get_data():
    return {};

def main():
    global epoch
    # If GPS is on, get Cartesian (position, velocity) vectors and UTC time.
    # Convert Cartesian coordinates and time to a Keplerian Elements array.
    # Update the config_adcs.yaml file with new KOE array.
    if(gps_is_on()):
        data = get_data() # Data is a list (cache) of dictionaries representing a timestep.
        i = len(data) # Get the last dictionary in the cache.
        r = [data[i]['x_pos'], data[i]['y_pos'], data[i]['z_pos']] # Position state vector.
        vel = [data[i]['x_vel'], data[i]['y_vel'], data[i]['z_vel']] # Velocity state vector.
        epoch = data[i]['time'] # Datetime object representing the epoch.
            
        drag = 2.2 # Fixed drag coefficient.
        koe_array = cart2kep(r, vel) # Convert state vectors into an array representing the KOE.
        koe_array = np.append(koe_array, drag) # 
        koe_array = np.insert(koe_array, 0, utc2jul(epoch))
        print(koe_array)
    # If GPS is off, write data to the YAML from the previous TLE file.
    # 
    else:
        write_config('config_adcs.yaml', tle_get_data())
        config = load_config('config_adcs.yaml')
        epoch = datetime.utcnow()
        koe_array = np.array([])
        for key, val in config['adcs']['koe'].items():
            koe_array = np.append(koe_array, val)
    # print(koe_array)
    # config['adcs']['sc']['jd0'] = utc2jul(epoch) # Use write config.

if __name__ == "__main__":
    t1 = threading.Thread(target=main, args=(), daemon=True)
    t1.start()
    t1.join()
    print("Calculation complete for " + epoch.strftime("%Y-%m-%d %H:%M:%S") + " UTC.")
