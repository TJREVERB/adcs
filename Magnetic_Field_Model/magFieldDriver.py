#Created by Bharath Dileepkumar
#Last Edited: December 6, 2018

import decyear, getDCM, sun_vec, utc2jul, wrldmagm, datetime
import numpy as np
import pymap3d as pm
class magFieldDriver:
     def testFunction(self, bV, sV, lat, long, alt, time):
         #epoch_date = datetime.datetime(1900, 1, 1, 0, 0)
         jdays = utc2jul.utc2jul(time) #- utc2jul.utc2jul(epoch_date)
         sI = sun_vec.sun_vec(jdays-2444244.5)
         print("sI:")
         print(sI)
         print("-----------------------------------------")
         print()
         gm = wrldmagm.wrldmagm("WMM2015.COF")
         bI = np.matrix(gm.wrldmagm(lat, long, alt, decyear.decyear(time)))
         print("bI:")
         print(bI)
         print("-----------------------------------------")
         print()
         bI = np.squeeze(np.asarray(bI))
         bI = pm.ned2ecef(bI[0], bI[1], bI[2], lat, long, alt)
         bI = pm.ecef2eci(bI, time)
         bI = np.matrix(bI).transpose()
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

# time = datetime.datetime(2018, 11, 8, 12, 0, 0)
# gm = wrldmagm.wrldmagm("WMM2015.COF")
# bI = np.matrix(gm.wrldmagm(5, -135, 396856, decyear.decyear(time)))
# print("bI:")
# print(bI)
