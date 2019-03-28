from pathlib import Path

from .get_q_ref
from .get_dcm import get_dcm
from .kep_to_cart import kep_to_cart
from .dec_year import dec_year
from .sun_vec import sun_vec
from .sun_sensors import sun_sensors
from .utc_to_jul import utc_to_jul
from .wrldmagm import WrldMagM
from .cart_to_kep import cart_to_kep
from .get_mc import get_mc
from .dcm_to_q import dcm_to_q
from core import load_config

from . import gps_dummy
from . import tle_dummy
from . import tle_points

import numpy as np
from datetime import datetime, date
from pymap3d import ecef2eci

import logging

logger = logging.Logger("ADCS")


def gps_is_on():
    return True

def start():
    global epoch
    config = load_config()  # Load the data from the YAML.
    # If GPS is on, get Cartesian (position, velocity) vectors and UTC time from the GPS.
    # Convert Cartesian coordinates and time to a Keplerian Elements array.
    # Generate a new TLE using the KOE.

    if gps_is_on():  # If we ask for GPS coordinates and the GPS responds:
        # Data is a list (cache) of dictionaries representing one timestep.
        data = gps_dummy.get_data()
        if gps_dummy.data_is_valid(data):  # If the data is valid:
            i = len(data)-1  # Get the last dictionary in the cache.
            # Position state vector.
            r = [data[i]['x_pos'], data[i]['y_pos'], data[i]['z_pos']]
            # Velocity state vector.
            vel = [data[i]['x_vel'], data[i]['y_vel'], data[i]['z_vel']]
            epoch = data[i]['time']  # Datetime object representing the epoch.

            # Convert state vectors into an array representing the KOE.
            koe_array = cart_to_kep(r, vel)
            koe_list = koe_array.tolist()
            # koe_array = np.insert(koe_array, 0, epoch)  # Add the datetime object epoch to the beginning.
            koe_list.insert(0, epoch)
            # koe_array = np.append(koe_array, data['adcs']['tledata']['bstardrag'])  # Append the B-star drag coefficient
            koe_list.append(config['adcs']['sc']['bstardrag'])
            temp_tle = tle_points.propagate(koe_list)  # Generate the new TLE.

            # print(koe_list)
            # print(temp_tle)

            # tjreverbtle = open(config['adcs']['tlefiles']['tjreverb'], "w")  # Open the main TJREVERB TLE for writing.
            # tjreverbtle.write(temp_tle)  # Write the new TLE to TJREVERB TLE.
            # tjreverbtle.close()  # Close the file.

            # Backup the TLE data.
            backuptle = open(config['adcs']['tlefiles']['backup'], "w")
            backuptle.write(temp_tle)
            backuptle.close()

            lla = tle_dummy.get_lla(epoch)  # Pull LLA data from TJREVERB TLE.
        else:  # If the GPS is on but the data is invalid:
            epoch = datetime.utcnow()  # Set current time to the system time.
            # Uses PyOrbital to propogate the TLE using epoch, which returns its LLA.
            lla = tle_dummy.get_lla(epoch)

    # If GPS is off, use the system time and TJREVERB TLE to propogate the current LLA.
    else:  # If we ask for GPS coordinates and the GPS not respond:
        epoch = datetime.utcnow()  # Set current time to the system time.
        # Uses PyOrbital to propogate the TLE using epoch, which returns its LLA.
        lla = tle_dummy.get_lla(epoch)

    # needed to incremement revnum
    # print(tle_dummy.get_xyz(epoch)['xyz_pos'])
    # print(tle_dummy.get_xyz(epoch)['xyz_vel'])
    # poskep = cart_to_kep(tle_dummy.get_xyz(epoch)['xyz_pos'], tle_dummy.get_xyz(epoch)['xyz_vel'])
    # print(poskep)
    # if (poskep[4]>0 and config['adcs']['tledata']['oldargp']<=0):
    #     with open("config_adcs.yaml") as f:
    #         list_doc = yaml.load(f)
    #     #print(type(list_doc))
    #     list_doc['adcs']['tledata']['revnum']=list_doc['adcs']['tledata']['revnum']+1
    #     with open("config_adcs.yaml", "w") as f:
    #         yaml.dump(list_doc, f, default_flow_style=False)
    # config = load_config('config_adcs.yaml')
    # config['adcs']['tledata']['oldargp'] = poskep[4]

    # write_config('config_adcs.yaml', utc_to_jul(epoch))  # config['adcs']['sc']['jd0'] = utc_to_jul(epoch)
    # Instantiates the WrldMagM object.
    gm = WrldMagM((Path(__file__).parent.resolve() /
                   config['adcs']['wrldmagm']))

    # Calculate the magnetic field vector in ECEF. Altitude is multiplied to convert meters to feet.
    magECEF = gm.wrldmagm(lla['lat'], lla['lon'], lla['alt'], date.today())

    magECEF = np.squeeze(np.asarray(magECEF))
    magECI = ecef2eci(magECEF, epoch)

    # Magnetic field in inertial frame, converts teslas to nanoteslas.
    bI = 1.0*(10e-09) * magECI
    bI = bI/np.linalg.norm(bI)
    bI = np.asmatrix(bI)
    bI = bI.getH()

    # Sun vector in intertial frame.
    sI = sun_vec(utc_to_jul(epoch)-utc_to_jul(datetime(1980, 1, 6, 0, 0, 0)))
    sI = sI/np.linalg.norm(sI)  # Normalize sI.

    print(bI)
    print(sI)

    # bV and sV data are taken from the onboard magnetometer and sun_sensors.

    # bV and sV data are taken from the onboard magnetometer and sunsensors.

    # dcm = get_dcm(bV, sV, bI, sI)
    # q = dcm_to_q(dcm)
    # qref = get_q_ref.get_q_ref(poskep)                           
    # qerr = get_q_err(q,qref)                              
    # thetaerr = get_theta_err(qerr)
    # mmax = [.1,.1,.1]
    # mtrans = np.matrix([[1,0,0],[0,1,0],[0,0,1]])
    # ctorque = np.matrix([0,0,0])
    # magdip = get_mc(ctorque.getH(),bV,mmax,mtrans)
    # ctprod = np.cross(magdip,bV)

    # isisimtq.k_imtq_start_actuation_dipole(imtq_axis_data(magdip[0], magdip[1], magdip[2]), 800)
