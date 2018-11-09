import logging
import os
import time
from subprocess import call
from threading import Thread

import numpy as np
import math
from math import floor, pi, sin, cos
from datetime import datetime

import pynmea2
import serial

from core import config
from . import aprs

logger = logging.getLogger("ADCS")

def send(msg):
    msg += "\n"
    ser.write(msg.encode("utf-8"))

"""""""""""""""""""""""
MAIN METHODS
"""""""""""""""""""""""
def utc2jul(dt):
    a = math.floor((14-dt.month)/12)
    y = dt.year + 4800 - a
    m = dt.month + 12*a - 3
    jdn = dt.day + math.floor((153*m + 2)/5) + 365*y + math.floor(y/4) - \
    math.floor(y/100) + math.floor(y/400) - 32045
    jd = jdn + (dt.hour - 12) / 24 + dt.minute / 1440 + \
    dt.second / 86400 - 2415020.5 #Subtracting dates since Jan 1, 1900
    # + dt.microsecond / 86400000000
    return jd


"""""""""""""""""""""""
DRIVING METHODS
"""""""""""""""""""""""
def listen():
    while True:
        epoch = np.matrix([2018, 11, 8, 12, 00, 00])
        time.sleep(1);

def updateVals(msg):
    #saves velocity data from gps
    global velocity_data
    velocity_data = msg

def keyin():
    while (True):
        # GET INPUT FROM YOUR OWN TERMINAL
        # TRY input("shihaoiscoolforcommentingstuff") IF raw_input() doesn't work
        in1 = input("Type Command: ")
        send(in1)
        # send("TJ" + in1 + chr(sum([ord(x) for x in "TJ" + in1]) % 128))


def on_startup():
    # GLOBAL VARIABLES ARE NEEDED IF YOU "CREATE" VARIABLES WITHIN THIS METHOD
    # AND ACCESS THEM ELSEWHERE
    global mainthread, logfile, tlt
    # cached_nmea_obj = (None,None)

    # serialPort = config['adcs']['serial_port']
    # REPLACE WITH COMx IF ON WINDOWS
    # REPLACE WITH /dev/ttyUSBx if 1 DOESNT WORK
    # serialPort = "/dev/ttyS3"
    # OPENS THE SERIAL PORT FOR ALL METHODS TO USE WITH 19200 BAUD
    # ser = serial.Serial(serialPort, 9600)
    # CREATES A THREAD THAT RUNS THE LISTEN METHOD
    mainthread = Thread(target=listen, args=(), daemon=True)
    #gpsdata = Thread(target=updateVals, args=(), daemon=True)
    mainthread.start()
    #gpsdata.start()

    tlt = time.localtime()

    # Open the log file
    log_dir = os.path.join(config['core']['log_dir'], 'adcs')
    filename = 'adcs' + '-'.join([str(x) for x in time.localtime()[0:3]])
    # ensure that the GPS log directory exists
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    logfile = open(os.path.join(log_dir, filename + '.txt'), 'a+')

    log('RUN@' + '-'.join([str(x) for x in tlt[3:5]]))

    # send("ANTENNAPOWER OFF")


# I NEED TO KNOW WHAT NEEDS TO BE DONE IN NORMAL, LOW POWER, AND EMERGENCY MODES
def enter_normal_mode():
    # UPDATE GPS MODULE INTERNAL COORDINATES EVERY 10 MINUTES
    # time.sleep(600)
    pass


def enter_low_power_mode():
    # UPDATE GPS MODULE INTERNAL COORDINATES EVERY HOUR
    # time.sleep(3600)
    pass


def enter_emergency_mode():
    pass

# TODO fix this
"""
def get_pry():
    return (-1,-1,-1)

def get_mag():
    return (-1,-1,-1)

def get_abs():
    return (-1,-1,-1)

def can_TJ_be_seen():
    return True # fix me!
"""
# USE THIS LOG FUNCTION
def log(msg):
    global logfile
    logfile.write(msg + '\n')
    logfile.flush()


if __name__ == '__main__':

    t2 = Thread(target=keyin, args=())
    t2.daemon = True
    t2.start()
    while True:
        time.sleep(1)
