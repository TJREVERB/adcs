# Orbital Decay
# Algorithm by Australian Government Bureau of Meteorology
# Translated from QuickBasic by Ayush Rautwar & Derek Goh
# Note: time increment should be below 5 days for the best accuracy

import math
from decimal import *
getcontext().prec = 100

def orbitalDecay():
   print("Enter the following satellite data:")
   n = input("Name:")
   m = float(input("Mass(kg):"))
   a = float(input("Area(square meters):"))
   h = float(input("Starting Height(km):"))
   f10 = float(input("Solar Radio Flux(SFU):"))
   ap = float(input("Geomagnetic A Index:"))
   
   print
   print
   print ("Satellite Name:", n)
   print
   print("Mass:", m,"kg")
   print("Area:", a,"square meters")
   print("Starting Height:", h,"km")
   print("Solar Radio Flux:",f10,"SFU" )
   print("Geomagnetic A Index:", ap)

   re, me, g = 6371000, (5.98*(10**24)), (6.67*(10**-11))
   t, dt = 0, 1 # t is time from launch, dt is time increment in days
   d9 = dt*3600*24
   h1=0
   h2=h
   r = (re+(h*1000))
   p = (2*math.pi*(((r**3)/me/g)**0.5))

   sh,dn,dp = 0,0,0
   
   print
   print  ("TIME     HEIGHT    Period    ")
   print  ("(Days)    (km)     (mins)   ")

 
   
   while h>=100:
      sh = (900+2.5*(f10-70)+1.5*ap)/(27-(0.012*(h-200)))
      dn = ((6*(10**-10))*(math.e**(-(h-175)/sh))) 
      dp = 3*math.pi*a/m*r*dn*d9
      pm = (p/60) 
      mm = (1440/pm)
      nmm = (1440/((p-dp))/60)
      decay = (dp/dt/p)*mm
          print('{:<10}'.format(t),'{:^10}'.format(h),'{:>10}'.format(p/60) )
      t += dt
      p -= dp
      r = (g*me*p*p/4/math.pi/math.pi)**(1/3)
      h = (r-re)/1000
      
      
   print ("Reentry after", t, "Days, or", t/365, "years")

orbitalDecay()

