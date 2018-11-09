import logging
import os
import time
from subprocess import call
from threading import Thread

import numpy as np
import math
from math import floor, pi, sin, cos, pi, radians
from datetime import datetime, utcnow, date
from orbital import earth, KeplerianElements, utilities
from orbital.utilities import Position, Velocity
from pymap3d import ecef2eci

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
    dt.second / 86400 - 2415020.5
    return jd
def sun_vec(start_day):
    jd = start_day
    L = (279.696678 + 0.9856473354*jd + 2.267e-13*(jd**2))
    Ms_r = (pi/180)*(358.475845 + 0.985600267*jd - (1.12e-13)*(jd**2) - \
    (7e-20)*(jd**3))
    dL = 1.918*sin(Ms_r) + 0.02*sin(2*Ms_r)
    L_sun = ((pi/180)*(L+dL))%(2*pi)
    inc_E = (pi/180)*(-23.45)
    #R = [1,0,0; 0,cos(inc_E),sin(inc_E); 0,-sin(inc_E),cos(inc_E)];# [3,3]
    #sun_ecl = [cos(L_sun);sin(L_sun);zeros(1,size(start_day,2))];  # [3,n]
    R = np.array([[1,0,0],[0,cos(inc_E),sin(inc_E)], \
    [0,-sin(inc_E),cos(inc_E)]],np.float32)
    sun_ecl = np.array([[cos(L_sun)],[sin(L_sun)],[0]],np.float32)
    sun_equ = np.matmul(R,sun_ecl)   # [3,1]
    return sun_equ
def kep2cart(KOE):
    orbitx = KeplerianElements(a=KOE[0], e=KOE[1], i=KOE[2], raan=KOE[3],
    arg_pe=KOE[4], M0=KOE[5], body=earth)
    cart = np.array([list(orbitx.r), list(orbitx.v)])
    return cart
def jd2dvec(jd):
    ps = jd - 2400000.5 # Done to increase time precision
    epochvec = list(jd2gcal(2400000.5, ps))  #Converts tuple to list
    hours = int(epochvec[3]*24)
    epochvec.append(epochvec[3]*24 - hours) #Sets jdarray[4] to decimal of hours
    epochvec[3] = hours
    minutes = int(epochvec[4]*60)
    epochvec.append(epochvec[4]*60 - minutes)
    epochvec[4] = minutes
    epochvec[5] = (epochvec[5]*60)
    return epochvec
def wrldmagm(self, dlat, dlon, h, time=date.today()):
    time = time.year+((time - date(time.year,1,1)).days/365.0)
    alt = h/3280.8399

    otime = oalt = olat = olon = -1000.0

    dt = time - self.epoch
    glat = dlat
    glon = dlon
    rlat = math.radians(glat)
    rlon = math.radians(glon)
    srlon = math.sin(rlon)
    srlat = math.sin(rlat)
    crlon = math.cos(rlon)
    crlat = math.cos(rlat)
    srlat2 = srlat*srlat
    crlat2 = crlat*crlat
    self.sp[1] = srlon
    self.cp[1] = crlon

    if (alt != oalt or glat != olat):
        q = math.sqrt(self.a2-self.c2*srlat2)
        q1 = alt*q
        q2 = ((q1+self.a2)/(q1+self.b2))*((q1+self.a2)/(q1+self.b2))
        ct = srlat/math.sqrt(q2*crlat2+srlat2)
        st = math.sqrt(1.0-(ct*ct))
        r2 = (alt*alt)+2.0*q1+(self.a4-self.c4*srlat2)/(q*q)
        r = math.sqrt(r2)
        d = math.sqrt(self.a2*crlat2+self.b2*srlat2)
        ca = (alt+d)/r
        sa = self.c2*crlat*srlat/(r*d)

    if (glon != olon):
        for m in range(2,self.maxord+1):
            self.sp[m] = self.sp[1]*self.cp[m-1]+self.cp[1]*self.sp[m-1]
            self.cp[m] = self.cp[1]*self.cp[m-1]-self.sp[1]*self.sp[m-1]

    aor = self.re/r
    ar = aor*aor
    br = bt = bp = bpp = 0.0
    for n in range(1,self.maxord+1):
        ar = ar*aor
        m=0
        D3=
        D4=(n+m+1)

        while D4>0:
            if (alt != oalt or glat != olat):
                if (n == m):
                    self.p[m][n] = st * self.p[m-1][n-1]
                    self.dp[m][n] = st*self.dp[m-1][n-1]+ct*self.p[m-1][n-1]

                elif (n == 1 and m == 0):
                    self.p[m][n] = ct*self.p[m][n-1]
                    self.dp[m][n] = ct*self.dp[m][n-1]-st*self.p[m][n-1]

                elif (n > 1 and n != m):
                    if (m > n-2):
                        self.p[m][n-2] = 0
                    if (m > n-2):
                        self.dp[m][n-2] = 0.0
                    self.p[m][n] = ct*self.p[m][n-1]-self.k[m][n]*self.p[m][n-2]
                    self.dp[m][n] = ct*self.dp[m][n-1] - st*self.p[m][n-1]-self.k[m][n]*self.dp[m][n-2]

            if (time != otime):
                self.tc[m][n] = self.c[m][n]+dt*self.cd[m][n]
                if (m != 0):
                    self.tc[n][m-1] = self.c[n][m-1]+dt*self.cd[n][m-1]
            par = ar*self.p[m][n]

            if (m == 0):
                temp1 = self.tc[m][n]*self.cp[m]
                temp2 = self.tc[m][n]*self.sp[m]
            else:
                temp1 = self.tc[m][n]*self.cp[m]+self.tc[n][m-1]*self.sp[m]
                temp2 = self.tc[m][n]*self.sp[m]-self.tc[n][m-1]*self.cp[m]

            bt = bt-ar*temp1*self.dp[m][n]
            bp = bp + (self.fm[m] * temp2 * par)
            br = br + (self.fn[n] * temp1 * par)

            if (st == 0.0 and m == 1):
                if (n == 1):
                    self.pp[n] = self.pp[n-1]
                else:
                    self.pp[n] = ct*self.pp[n-1]-self.k[m][n]*self.pp[n-2]
                parp = ar*self.pp[n]
                bpp = bpp + (self.fm[m]*temp2*parp)
            D4=D4-1
            m=m+1
    if (st == 0.0):
        bp = bpp
    else:
        bp = bp/st

    bx = -bt*ca-br*sa
    by = bp
    bz = bt*sa-br*ca

    bh = math.sqrt((bx*bx)+(by*by))
    ti = math.sqrt((bh*bh)+(bz*bz))
    dec = math.degrees(math.atan2(by,bx))
    dip = math.degrees(math.atan2(bz,bh))

    gv = -999.0
    if (math.fabs(glat) >= 55.):
        if (glat > 0.0 and glon >= 0.0):
            gv = dec-glon
        if (glat > 0.0 and glon < 0.0):
            gv = dec+math.fabs(glon);
        if (glat < 0.0 and glon >= 0.0):
            gv = dec+glon
        if (glat < 0.0 and glon < 0.0):
            gv = dec-math.fabs(glon)
        if (gv > +180.0):
            gv = gv - 360.0
        if (gv < -180.0):
            gv = gv + 360.0

    otime = time
    oalt = alt
    olat = glat
    olon = glon

"""""""""""""""""""""""
DRIVING METHODS
"""""""""""""""""""""""
def listen():
    while True:
        epoch = datetime.utcnow()
        config['adcs']['sc']['jd0'] = utc2jul(epoch)
        config['adcs']['sc']['inertia'] = np.diag([.0108, .0108, .0108])
        #Orbital Properties
        GM = 3.986004418*(10**14)
        KOE = np.array(config['adcs']['KOE']['sma'], \
        config['adcs']['KOE']['ecc'], config['adcs']['KOE']['incl'] \
        config['adcs']['KOE']['argp'], config['adcs']['KOE']['tran'])
        current_time = epoch
        epoch = np.array(epoch.year, epoch.month, epoch.day, epoch.hour, \
        epoch.minute, epoch.second)
        cart = kep2cart(KOE)
        cartloc = np.array(cart[0], cart[1], cart[2])
        #Magnetic Field Model
        epochvec = jd2dvec(config['adcs']['sc']['jd0'])
        lla = np.array(gps.lon gps.lat, gps.alt)
        magECEF = wrldmagm(lla(3),lla(1),lla(2), \
        decyear('01-January-2018','dd-mm-yyyy'))
        magECI = ecef2eci(magECEF, current_time)
        time.sleep(1);

def updateVals(msg):
    #saves velocity data from gps
    #If status is not 0, something is wrong
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
