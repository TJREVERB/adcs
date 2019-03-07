# Main ADCS Driver
# Python Methods Used:

from getDCM import getDCM
from kep2cart import kep2cart
from decyear import decyear
from kepel import kepel
from q2dcm import q2dcm
from sun_vec import sun_vec
from sunsensors import sunsensors
from utc2jul import utc2jul
from wrldmagm import wrldmagm
from cart2kep import cart2kep

import gps_dummy
import tle_dummy
import tle_points

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
    return False

def tle_get_data():
    return {}

def generate_tle(koe):
    return {}

def main():
    global epoch
    config = load_config('config_adcs.yaml')  # Load the data from the YAML.
    # If GPS is on, get Cartesian (position, velocity) vectors and UTC time from the GPS.
    # Convert Cartesian coordinates and time to a Keplerian Elements array.
    # Generate a new TLE using the KOE.

    if gps_is_on():  # If we ask for GPS coordinates and the GPS responds:
        data = gps_dummy.get_data()  # Data is a list (cache) of dictionaries representing one timestep.
        if gps_dummy.data_is_valid(data):  # If the data is valid:
            i = len(data)-1  # Get the last dictionary in the cache.
            r = [data[i]['x_pos'], data[i]['y_pos'], data[i]['z_pos']]  # Position state vector.
            vel = [data[i]['x_vel'], data[i]['y_vel'], data[i]['z_vel']]  # Velocity state vector.
            epoch = data[i]['time']  # Datetime object representing the epoch.

            koe_array = cart2kep(r, vel)  # Convert state vectors into an array representing the KOE.
            koe_array = np.insert(koe_array, 0, epoch)  # Add the datetime object epoch to the beginning.
            koe_array = np.append(koe_array, data['adcs']['tledata']['bstardrag'])  # Append the B-star drag coefficient
            temp_tle = tle_points.propagate(koe_array)  # Generate the new TLE.

            tjreverbtle = open("tjreverb_tle.txt", "w")  # Open the main TJREVERB TLE for writing.
            tjreverbtle.write(temp_tle)  # Write the new TLE to TJREVERB TLE.
            tjreverbtle.close()  # Close the file.

            backuptle = open("backup_tle.txt", "w")  # Backup the TLE data.
            backuptle.write(temp_tle)
            backuptle.close()

            lla = tle_dummy.get_lla(epoch)  # Pull LLA data from TJREVERB TLE.
        else:  # If the GPS is on but the data is invalid:
            epoch = datetime.utcnow()  # Set current time to the system time.
            lla = tle_dummy.get_lla(epoch)  # Uses PyOrbital to propogate the TLE using epoch, which returns its LLA.

    # If GPS is off, use the system time and TJREVERB TLE to propogate the current LLA.
    else:  # If we ask for GPS coordinates and the GPS not respond:
        epoch = datetime.utcnow()  # Set current time to the system time.
        lla = tle_dummy.get_lla(epoch)  # Uses PyOrbital to propogate the TLE using epoch, which returns its LLA.
        # Run Ayush's code to get the propogated TLE.

        # write_config('config_adcs.yaml', tle_get_data()) # Write TLE data to YAML.

        """
        epoch = datetime.utcnow()
        koe_array = np.array([utc2jul(epoch)])
        for key, val in config['adcs']['koe'].items():
            koe_array = np.append(koe_array, val)
        """

    # write_config('config_adcs.yaml', utc2jul(epoch))  # config['adcs']['sc']['jd0'] = utc2jul(epoch)
    gm = wrldmagm(config['adcs']['wrldmagm'])  # Instantiates the wrldmagm object.
    magECEF = gm.wrldmagm(lat, lon, alt, decyear(datetime(2018, 1, 1)))


if __name__ == "__main__":
    t1 = threading.Thread(target=main, args=(), daemon=True)
    t1.start()
    t1.join()
    print("Calculation complete for " + epoch.strftime("%Y-%m-%d %H:%M:%S") + " UTC.")
