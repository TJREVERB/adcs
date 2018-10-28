import numpy as np
def cubesatSim():

    #Spacecraft Properties
    epoch = np.matrix([2018, 1, 1, 12, 00, 00])
    sc.jd0=juliandate(epoch)
    #sc.mjd0 = mjuliandate(epoch)
    sc.inertia=np.diag([.0108, .0108, .0108])
    #sc.rc = 6771                                      
    #sc.In = 51.6*pi/180                                
    #c.Torb = 2*pi/sqrt(3.986e5/sc.rc^3)

    #Orbital Properties
    GM = 3.986004418*(10**14)
    KOE.sma = 6778557                                
    KOE.ecc = 0.0001973                                
    KOE.incl = 51.6397*pi/180                       
    KOE.argp = 211.4532*pi/180                          
    KOE.tran = 0*pi/180
    koeVect = struct2array(KOE)                       
    epoch = [2018, 1, 1, 12, 00, 00]                 
    cart = kep2cart(KOE)                              
    for i in range(0,2):
        cartloc[i] = cart[i]                               
    #lla = eci2lla(cartloc.getH(),epoch)

    #Magnetic Field Model
    epochvec = datevec(datetime(sc.jd0,'convertfrom','juliandate')) #Jason, please do the magic
    lla = eci2lla(cartloc.getH(),epochvec)
    eci2ecef = dcmeci2ecef('IAU-2000/2006',epochvec) #needs to be written
    ecef2eci = eci2ecef.getH()
    magECEF = wrldmagm(lla(3),lla(1),lla(2),decyear('01-January-2018','dd-mm-yyyy'))
    magECI = ecef2eci*magECEF 
    
