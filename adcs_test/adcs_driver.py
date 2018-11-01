import logging
import os
import time
from subprocess import call
from threading import Thread

import pynmea2
import serial

from core import config
from . import aprs
from submodules import gps
from adcs import sun_vec

logger = logging.getLogger("ADCS")

def send(msg):
    msg += "\n"
    ser.write(msg.encode("utf-8"))


def listen():
    #There may be delays into the computations, try parallel threads
    #TEST WITH PROCESSOR
    while True:
        sun_vec.sun_vec(time.time()/60/60/24)
        time.sleep(1)

def on_startup():
    # GLOBAL VARIABLES ARE NEEDED IF YOU "CREATE" VARIABLES WITHIN THIS METHOD
    # AND ACCESS THEM ELSEWHERE
    global dcmThread, logfile
    # cached_nmea_obj = (None,None)

    # serialPort = config['adcs']['serial_port']
    # REPLACE WITH COMx IF ON WINDOWS
    # REPLACE WITH /dev/ttyUSBx if 1 DOESNT WORK
    # serialPort = "/dev/ttyS3"
    # OPENS THE SERIAL PORT FOR ALL METHODS TO USE WITH 19200 BAUD
    # ser = serial.Serial(serialPort, 9600)
    # CREATES A THREAD THAT RUNS THE LISTEN METHOD
    dcmThread = Thread(target=listen, args=(), daemon=True)
    dcmThread.start()

    # Open the log file
    log_dir = os.path.join(config['core']['log_dir'], 'adcs')
    filename = 'adcs' + '-'.join([str(x) for x in time.localtime()[0:3]])
    # ensure that the GPS log directory exists
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    logfile = open(os.path.join(log_dir, filename + '.txt'), 'a+')

    log('RUN@' + '-'.join([str(x) for x in tlt[3:5]]))

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
