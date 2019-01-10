#Created by Bharath Dileepkumar
#Last Edited: January 10, 2018

import decyear, getDCM, sun_vec, utc2jul, wrldmagm, datetime
import numpy as np
import pymap3d as pm
import math
class magFieldDriver:
     def testFunction(self, bV, sV, lat, long, alt, time):
         jdays = utc2jul.utc2jul(time)
         sI = sun_vec.sun_vec(jdays-2444244.5)
         temp = np.matrix([0, 0, -1]).transpose()
         temp = pm.ned2ecef(0, 0, -1, 45, 30, 0)
         temp = temp/np.linalg.norm(temp)
         print(temp)
         print("sI:")
         print(sI)
         print("-----------------------------------------")
         print()
         gm = wrldmagm.wrldmagm("WMM2015.COF")
         bI = np.matrix(gm.wrldmagm(lat, long, alt, decyear.decyear(time)))
         bI = bI/np.linalg.norm(bI)
         print("bI: (NED)")
         print(bI)
         print("-----------------------------------------")
         print()
         bI = np.squeeze(np.asarray(bI))
         bI = pm.ned2ecef(bI[0], bI[1], bI[2], lat, long, alt)
         bIECEF = np.matrix(bI).transpose()
         bIECEF = bI/np.linalg.norm(bI)
         print("bI: (ECEF)")
         print(bIECEF)
         print("-----------------------------------------")
         print()
         bI = pm.ecef2eci(bI, time)
         bI = np.matrix(bI).transpose()
         bI = bI/np.linalg.norm(bI)
         print("bI: (ECI)")
         print(bI)
         print("------------------------------")
         print()
         print("DCM:")
         return getDCM.getDCM(bV, sV, bI, sI)

     def mainTester(self):
         new_date = datetime.datetime(2018, 11, 8, 12, 0, 0)
         sV = np.matrix([0.736594581345171, -0.321367778737346, 0.595106018724694]).transpose()
         bV = np.matrix([ 0.593302194154829, -0.297353472232347, -0.748046401610510]).transpose()
         lat_deg = 5.335745187657780
         lon_deg = -1.348386750055788e+02
         alt_meter = 3.968562753276266e+05
         return self.testFunction(bV, sV, lat_deg, lon_deg, alt_meter, new_date)

driver = magFieldDriver()
print(driver.mainTester())

