# geomag.py
# by Christopher Weiss cmweiss@gmail.com
# Modified by Bharath Dileepkumar, 11-08-2018

# Adapted from the geomagc software and World Magnetic Model of the NOAA
# Satellite and Information Service, National Geophysical Data Center
# http://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml
#
# Suggestions for improvements are appreciated.

# USAGE:
#
# >>> gm = geomag.GeoMag("WMM.COF")
# >>> mag = gm.GeoMag(80,0)
# >>> mag.dec
# -6.1335150785195536
# >>>

import math, os, unittest
import numpy as np
import pymap3d as pm
from datetime import date

class wrldmagm:

    def wrldmagm(self, dlat, dlon, h, time): # latitude (decimal degrees), longitude (decimal degrees), altitude (feet), date
        #time = date('Y') + date('z')/365
        #time = time.year+((time - date(time.year,1,1)).days/365.0)
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

        #/* CONVERT FROM GEODETIC COORDS. TO SPHERICAL COORDS. */
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

            #for (m=0,D3=1,D4=(n+m+D3)/D3;D4>0;D4--,m+=D3):
            m=0
            D3=1
            #D4=(n+m+D3)/D3
            D4=(n+m+1)
            while D4>0:

        # /*
                # COMPUTE UNNORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
                # AND DERIVATIVES VIA RECURSION RELATIONS
        # */
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

        # /*
                # TIME ADJUST THE GAUSS COEFFICIENTS
        # */
                if (time != otime):
                    self.tc[m][n] = self.c[m][n]+dt*self.cd[m][n]
                    if (m != 0):
                        self.tc[n][m-1] = self.c[n][m-1]+dt*self.cd[n][m-1]

        # /*
                # ACCUMULATE TERMS OF THE SPHERICAL HARMONIC EXPANSIONS
        # */
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
        # /*
                    # SPECIAL CASE:  NORTH/SOUTH GEOGRAPHIC POLES
        # */
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
        # /*
            # ROTATE MAGNETIC VECTOR COMPONENTS FROM SPHERICAL TO
            # GEODETIC COORDINATES
        # */
        bx = -bt*ca-br*sa
        by = bp
        bz = bt*sa-br*ca
        # /*
            # COMPUTE DECLINATION (DEC), INCLINATION (DIP) AND
            # TOTAL INTENSITY (TI)
        # */
        bh = math.sqrt((bx*bx)+(by*by))
        ti = math.sqrt((bh*bh)+(bz*bz))
        dec = math.degrees(math.atan2(by,bx))
        dip = math.degrees(math.atan2(bz,bh))
        # /*
            # COMPUTE MAGNETIC GRID VARIATION IF THE CURRENT
            # GEODETIC POSITION IS IN THE ARCTIC OR ANTARCTIC
            # (I.E. GLAT > +55 DEGREES OR GLAT < -55 DEGREES)

            # OTHERWISE, SET MAGNETIC GRID VARIATION TO -999.0
        # */
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

        class RetObj:
            pass
        retobj = RetObj()
        retobj.dec = dec
        retobj.dip = dip
        retobj.ti = ti
        retobj.bh = bh
        retobj.bx = bx
        retobj.by = by
        retobj.bz = bz
        retobj.lat = dlat
        retobj.lon = dlon
        retobj.alt = h
        retobj.time = time
        retMag = np.matrix([retobj.bx, retobj.by, retobj.bz])
        retMag = retMag.transpose()
        #retFinal = np.matrix([retMag, retobj.bh, dec, retobj.dip, retobj.ti])
        #return retobj
        return retMag

    def __init__(self, wmm_filename=None):
        if not wmm_filename:
            wmm_filename = os.path.join(os.path.dirname(__file__), 'WMM.COF')
        wmm=[]
        with open(wmm_filename) as wmm_file:
            for line in wmm_file:
                linevals = line.strip().split()
                if len(linevals) == 3:
                    self.epoch = float(linevals[0])
                    self.model = linevals[1]
                    self.modeldate = linevals[2]
                elif len(linevals) == 6:
                    linedict = {'n': int(float(linevals[0])),
                    'm': int(float(linevals[1])),
                    'gnm': float(linevals[2]),
                    'hnm': float(linevals[3]),
                    'dgnm': float(linevals[4]),
                    'dhnm': float(linevals[5])}
                    wmm.append(linedict)

        z = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        self.maxord = self.maxdeg = 12
        self.tc = [z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13]]
        self.sp = z[0:14]
        self.cp = z[0:14]
        self.cp[0] = 1.0
        self.pp = z[0:13]
        self.pp[0] = 1.0
        self.p = [z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14]]
        self.p[0][0] = 1.0
        self.dp = [z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13]]
        self.a = 6378.137
        self.b = 6356.7523142
        self.re = 6371.2
        self.a2 = self.a*self.a
        self.b2 = self.b*self.b
        self.c2 = self.a2-self.b2
        self.a4 = self.a2*self.a2
        self.b4 = self.b2*self.b2
        self.c4 = self.a4 - self.b4

        self.c = [z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14]]
        self.cd = [z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14]]

        for wmmnm in wmm:
            m = wmmnm['m']
            n = wmmnm['n']
            gnm = wmmnm['gnm']
            hnm = wmmnm['hnm']
            dgnm = wmmnm['dgnm']
            dhnm = wmmnm['dhnm']
            if (m <= n):
                self.c[m][n] = gnm
                self.cd[m][n] = dgnm
                if (m != 0):
                    self.c[n][m-1] = hnm
                    self.cd[n][m-1] = dhnm

        #/* CONVERT SCHMIDT NORMALIZED GAUSS COEFFICIENTS TO UNNORMALIZED */
        self.snorm = [z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13]]
        self.snorm[0][0] = 1.0
        self.k = [z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13]]
        self.k[1][1] = 0.0
        self.fn = [0.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0]
        self.fm = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0]
        for n in range(1,self.maxord+1):
            self.snorm[0][n] = self.snorm[0][n-1]*(2.0*n-1)/n
            j=2.0
            #for (m=0,D1=1,D2=(n-m+D1)/D1;D2>0;D2--,m+=D1):
            m=0
            D1=1
            D2=(n-m+D1)/D1
            while (D2 > 0):
                self.k[m][n] = (((n-1)*(n-1))-(m*m))/((2.0*n-1)*(2.0*n-3.0))
                if (m > 0):
                    flnmj = ((n-m+1.0)*j)/(n+m)
                    self.snorm[m][n] = self.snorm[m-1][n]*math.sqrt(flnmj)
                    j = 1.0
                    self.c[n][m-1] = self.snorm[m][n]*self.c[n][m-1]
                    self.cd[n][m-1] = self.snorm[m][n]*self.cd[n][m-1]
                self.c[m][n] = self.snorm[m][n]*self.c[m][n]
                self.cd[m][n] = self.snorm[m][n]*self.cd[m][n]
                D2=D2-1
                m=m+D1

class GeoMagTest(unittest.TestCase):

    d1=date(2015,1,1)
    d2=date(2017,7,2)

    test_values = (
        # date, alt, lat, lon, var
        (d1, 0, 80, 0,  -3.85),
        (d1, 0, 0, 120, 0.57),
        (d1, 0, -80, 240,  69.81),
        (d1, 328083.99, 80, 0, -4.27),
        (d1, 328083.99, 0, 120, 0.56),
        (d1, 328083.99, -80, 240, 69.22),
        (d2, 0, 80, 0, -2.75),
        (d2, 0, 0, 120, 0.32),
        (d2, 0, -80, 240, 69.58),
        (d2, 328083.99, 80, 0, -3.17),
        (d2, 328083.99, 0, 120, 0.32),
        (d2, 328083.99, -80, 240, 69.00),
    )

    # def test_declination(self):
    #     gm = GeoMag()
    #     for values in self.test_values:
    #         calcval=gm.GeoMag(values[2], values[3], values[1], values[0])
    #         print(calcval)
    #         self.assertAlmostEqual(values[4], calcval.dec, 2, 'Expected %s, result %s' % (values[4], calcval.dec))

#if __name__ == '__main__':
#    unittest.main()

def q2dcm(q):
    R = np.zeros((3,3))

    R[0,0] = q[0]^2-q[1]^2-q[2]^2+q[3]^2
    R[0,1] = 2*(q[0]*q[1]+q[2]*q[3])
    R[0,2] = 2*(q[0]*q[2]-q[1]*q[3])

    R[1,0] = 2*(q[0]*q[1]-q[2]*q[3])
    R[1,1] = -q[0]^2+q[1]^2-q[2]^2+q[3]^2
    R[1,2] = 2*(q[1]*q[2]+q[0]*q[3])

    R[2,0] = 2*(q[0]*q[2]+q[1]*q[3])
    R[2,1] = 2*(q[1]*q[2]-q[0]*q[3])
    R[2,2] = -q[0]^2-q[1]^2+q[2]^2+q[3]^2

    return R

import math
from math import pi
from math import sin
from math import cos
def sun_vec(start_day):
    # Julian days since Jan 0,1900
    #  Reference for this calculation is JD 2,415,020 which
    #  corresponds to 12:00:00 Jan 0,1900 ET (or 12:00:00 Dec 31,1899)
    jd = 29224.5 + start_day
    #  Mean longitude of sun, measured in the ecliptic from mean
    #  equinox of date:
    L = (279.696678 + 0.9856473354*jd + 2.267e-13*(jd**2))
    #  Mean anomaly of sun in radians
    Ms_r = (pi/180)*(358.475845 + 0.985600267*jd - (1.12e-13)*(jd**2) - (7e-20)*(jd**3))
    #  Correction between mean longitude and true longitude
    dL = 1.918*sin(Ms_r) + 0.02*sin(2*Ms_r)
    #  True longitude of sun, in radians
    L_sun = ((pi/180)*(L+dL))%(2*pi)
    #  Compute sun unit vector in ECI frame, where the Earth's
    #  equatorial plane is inclined inc_E radians to the ecliptic
    #  R defines a rotation about the x-axis
    inc_E = (pi/180)*(-23.45)
    #R = [1,0,0; 0,cos(inc_E),sin(inc_E); 0,-sin(inc_E),cos(inc_E)];# [3,3] #CONVERT
    #sun_ecl = [cos(L_sun);sin(L_sun);zeros(1,size(start_day,2))];  # [3,n]

    R = np.array([[1,0,0],[0,cos(inc_E),sin(inc_E)],[0,-sin(inc_E),cos(inc_E)]],np.float32)
    sun_ecl = np.array([[cos(L_sun)],[sin(L_sun)],[0]],np.float32)
    #  Since R is constant through time, can do a simple matrix multiply
    sun_equ = np.matmul(R,sun_ecl)   # [3,1]
    return sun_equ

from math import floor
"""
https://github.com/dannyzed/julian/blob/master/julian/julian.py
"""
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
"""
a = math.floor((14-epoch[1])/12)
y = epoch[0] + 4800 - a
m = epoch[1] + 12*a - 3
jdn = epoch[2] + math.floor((153*m + 2)/5) + 365*y + math.floor(y/4) - \
math.floor(y/100) + math.floor(y/400) - 32045
jd = jdn + (epoch[3] - 12) / 24 + epoch[4] / 1440 + \
epoch[5] / 86400 - 2415020.5 #Subtracting dates since Jan 1, 1900
# + dt.microsecond / 86400000000
"""

def getDCM(bV, sV, bI, sI):
  bV = np.matrix(bV)
  sV = np.matrix(sV)
  bI = np.matrix(bI)
  sI = np.matrix(sI)

  bV = np.reshape(bV, (1,-1))/LA.norm(bV) #
  sV = np.reshape(sV, (1,-1))/LA.norm(sV)  #
  bI = np.reshape(bI, (1,-1))/LA.norm(bI)  #
  sI = np.reshape(sI, (1,-1))/LA.norm(sI)  #

  vu2 = np.asmatrix(np.cross(bV, sV))
  vu2 = np.asmatrix(vu2/LA.norm(vu2))
  vmV = np.hstack((bV.getH(), vu2.getH(), np.asmatrix(np.cross(bV, vu2)).getH())) #
  iu2 = np.asmatrix(np.cross(bI, sI))
  iu2 = np.asmatrix(iu2/LA.norm(iu2))
  imV = np.hstack((bI.getH(), iu2.getH(), np.asmatrix(np.cross(bI, iu2)).getH())) #
  ivDCM = np.asmatrix(vmV)*np.asmatrix(imV).getH()
  return ivDCM

def decyear(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

def testFunction(bV, sV, lat, long, alt, time):
    jdays = utc2jul(time) - utc2jul(datetime.datetime(1900, 1, 1))
    sI = sun_vec(jdays)
    gm = wrldmagm("WMM.COF")
    bI = gm.wrldmagm(5.335745187657780, -1.348386750055788e+02, 3.968562753276266e+05, 	2018.8534)#temporary)
    bV = pm.ecef2eci(bI)
    return getDCM(bV, sV, bI, sI)

import datetime
def mainTester():
    new_date = datetime.datetime(2018, 11, 8, 12, 0)
    sV = np.matrix([0.736594581345171, -0.321367778737346, 0.595106018724694]).transpose()
    bV = np.matrix([ 0.593302194154829, -0.297353472232347, -0.748046401610510]).transpose()
    lat_deg = 5.335745187657780
    lon_deg = -1.348386750055788e+02
    alt_meter = 3.968562753276266e+05
    return testFunction(bV, sV, lat_deg, lon_deg, alt_meter, new_date)

    
    
    
