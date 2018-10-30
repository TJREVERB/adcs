import numpy as np
import math
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

    #Initial CubeSat Attitude
    #qtrue = [.5,.5,.5,.99];
    #qtrue = [0.5435   -0.0028   -0.6124   -0.5741];
    qtrue = np.matrix([0,0,sqrt(2)/2,sqrt(2)/2])
    qtrue = qtrue/normalize(qtrue)
    DCMtrue = q2dcm(qtrue)

    #Sensor Outputs
    #[magTotal,~] = BDipole(cart,sc.jd0,[0;0;0]);
    bI = 1.0*(10**(-09)) * magECI
    bI = bI/normalize(bI)
    sI = sun_vec(sc.jd0-juliandate([1980,1,6]))
    sI = sI/normalize(sI)
    bV = DCMtrue*bI
    bV = bV/norm(bV)
    sV = DCMtrue*sI
    sV = sV/norm(sV)

    #Attitude properties
    dcm = getDCM(bV,sV,bI,sI)
    q = dcm2q(dcm)

    #qref = getqref(koeVect)
    #qerr = getqerr(q,qref)
    #thetaerr = getthetaerr(qerr)


    #Simulation parameters
    Ts = 100
    tcamp = 1
    tint = .1
    tt = 0:tcamp:Ts

    #Initial angular rates
    n=sqrt(GM/(KOE.sma^3))
    %w = [0,-n,0]
    w = [0,0,0]
    Yt = [q, w;zeros(length(tt)-1,7)]

    #Control law
    sim.gain = 2e-5
    sim.rgain = 2e-5
    sim.mmax = np.matrix([.2,.2,.2])
    sim.mtrans = np.matrix([1,0,0;0,1,0;0,0,1])
    sim.wref = w

    #magdip = getMC(ctcomm,bV,mmax,mtrans)
    #ctprod = cross(magdip,bV)


    #Inputs for plotting
    magmoment = zeros(3,length(tt))
    magfield = zeros(3,length(tt))
    atterr = zeros(3,length(tt))

    #Simulation

    for i in range (1,length(tt)-1):
        KOEt = kepprop2b(KOE,tt(i),GM)
        jd = sc.jd0 + (i*tcamp/86400)
        #sc.mjd0 = sc.mjd0 + (tcamp/86400);
        #Ys = ode4(@(t,Y)AttDyn(t,Y,sc,sim,KOEt),[tt(i):tint:tt(i)+tcamp],Yt(i,:).getH());
        [~,Ys] = ode4(@(t,Y)AttDyn(t,Y,sc,sim,KOEt,jd),[tt(i):tint:tt(i)+tcamp],Yt(i,:).getH())
        Yt(i+1,:)=Ys(end,:)
        Yt(i+1,1:4) = Yt(i+1,1:4)/norm(Yt(i+1,1:4))

    [~,magfield(:,1),magmoment(:,1),atterr(:,1)] = AttDyn(tt(1),Yt(1,:).getH(),sc,sim,KOE,sc.jd0);

    for j in range(2,length(tt)):
        KOEt = kepprop2b(KOE,tt(j-1),GM)
        jd = sc.jd0 + ((j-1)*tcamp/86400)
        [~,magfield(:,j),magmoment(:,j),atterr(:,j)] = AttDyn(tt(j),Yt(j,:).getH(),sc,sim,KOEt,jd)
        # Transform inertial attitude to orbital attitude
        # TBI = TBO*TOI ==> TBI*(TOI)^-1 = TBO
        # If we have inertial position and velocity from GPS, we can use a different formula
        #
        #q = qmult(Yt(j,1:4)',qinv(qmult(getq(1,-pi/2),qmult(getq(3,sc.wo*tt(i)+pi/2),getq(1,sc.inertia)))));
        #att(:,i) = 2*q(1:3);
