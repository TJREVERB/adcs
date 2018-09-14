%% Created by Anonto Zaman 1/18/2018
% cubesatSim
% Basic 2U CubeSat motion simulation
% Assumes no disturbance torques, no gravity gradient 

%% Spacecraft properties
epoch = [2018, 1, 1, 12, 00, 00];                    % Time in UTC (3/17/2018 12:00:00 UTC)mjd = mjuliandate([y,mo,d,h,mi,s])
sc.jd0 = juliandate(epoch);                          % Initial satellite Julian Date
%sc.mjd0 = mjuliandate(epoch);                        % Iniitial Modified Julian Date
sc.inertia = diag([.0108 .0108 .00426]);             % Assuming mass of 2.6kg and side length of .1x.1x.2 meters 
% sc.rc = 6771;                                        % Orbit radius (km) (LEO), based on ISS average altitude
%sc.In = 51.6*pi/180;                                 % Orbit inclination (rad), based on ISS inclination
%sc.Torb = 2*pi/sqrt(3.986e5/sc.rc^3);                % Orbital period (s), based on calculated orbital period from ISS altitude, note decay

%% Orbital properties
GM = 3.986004418e14;                                 % (m^3/s^2)
KOE.sma = 6778557;                                % Semi-major axis (m)
KOE.ecc = 0.0001973;                                 % Eccentricity
KOE.incl = 51.6397*pi/180;                           % Inclination (rad)
KOE.raan = 115.7654*pi/180;                          % Right ascension of the ascending node (rad)
KOE.argp = 211.4532*pi/180;                          % Argument of perigee (rad)
KOE.tran = 0*pi/180;                                 % True anomaly (rad)
koeVect = struct2array(KOE);                         % Keplerian orbital elements in vector form
epoch = [2018, 1, 1, 12, 00, 00];                    % Time in UTC (3/17/2018 12:00:00 UTC)
cart = kep2cart(KOE);                                % Cartesian coordinates of satellite [pos,vel]
cartloc = cart(1:3);                                 % Cartesian location of satellite (row vector form)
%lla = eci2lla(cartloc',epoch);                       % Converts cartesian (ECI) coordinates to Latitude, Longitude, Altitude coordinates

%% Magnetic Field Model
epochvec = datevec(datetime(sc.jd0,'convertfrom','juliandate'));
lla = eci2lla(cartloc',epochvec); 
eci2ecef = dcmeci2ecef('IAU-2000/2006',epochvec);
ecef2eci = eci2ecef';
magECEF = wrldmagm(lla(3),lla(1),lla(2),decyear('01-January-2018','dd-mm-yyyy'));
magECI = ecef2eci*magECEF;  

%% Initial CubeSat Attitude
%qtrue = [.5,.5,.5,.99];                     % Initial rotation quaternion (rotation from Inertial to Vehicle Frame)
%qtrue = [0.5435   -0.0028   -0.6124   -0.5741];
qtrue = [0,0,sqrt(2)/2,sqrt(2)/2];
qtrue = qtrue/norm(qtrue);         % Formats rotation quaternion
DCMtrue = q2dcm(qtrue);                  % Initial true quaternion in DCM form

%% Sensor Outputs
%[magTotal,~] = BDipole(cart,sc.jd0,[0;0;0]);        % Magnetic field in nT (inputs are location in ECI frame, Julian Date, and Initial Angular rate)
bI = 1.0e-09 * magECI;                               % Magnetic field vector in ECI frame (measured in T)
bI = bI/norm(bI);
sI = sun_vec(sc.jd0-juliandate([1980,1,6]));         % Sun vector in ECI frame
sI = sI/norm(sI);
bV = DCMtrue*bI;                               % Magnetic field vector in vehicle frame
bV = bV/norm(bV);
sV = DCMtrue*sI;                               % Sun vector in vehicle frame
sV = sV/norm(sV);

%% Attitude properties
dcm = getDCM(bV,sV,bI,sI);                           % DCM from vehicle to inertial frame
q = dcm2q(dcm);                                      % Spacecraft initial quaternion
%{
qref = getqref(koeVect);                             % Calculates reference quaternion
qerr = getqerr(q,qref);                              % Calcualtes error quaternion
thetaerr = getthetaerr(qerr);                        % Error angles for each axis 
%}

%% Simulation parameters
Ts = 100;                                           % Simulation time (sec)
tcamp = 1;                                          % Plot sampling time/integration step
tint = .1;                                            % Integration step 
tt = 0:tcamp:Ts;                                     % Time vector for simulation

%% Initial angular rates
n=sqrt(GM/(KOE.sma^3));
%w = [0,-n,0];
w = [0,0,0];
Yt = [q, w;zeros(length(tt)-1,7)];

%% Control law 
sim.gain = 2e-5;                                              % Test gain value
sim.rgain = 2e-5;
sim.mmax = [.2,.2,.2];                                      % Maximum magnetic moment for each torquer axis
sim.mtrans = [1,0,0;0,1,0;0,0,1];                        % DCM to convert from body frame to magnetorquer frame
sim.wref = w;
%{
magdip = getMC(ctcomm,bV,mmax,mtrans);               % Magnetic dipole in body frame
ctprod = cross(magdip,bV);                           % Produced control torque
%}

%% Inputs for plotting
magmoment = zeros(3,length(tt));
magfield = zeros(3,length(tt));
atterr = zeros(3,length(tt));

%% Simulation

for i = 1:length(tt)-1
    KOEt = kepprop2b(KOE,tt(i),GM);                         % Propogates orbit
    jd = sc.jd0 + (i*tcamp/86400);                        % Updates epoch (Julian Date)
    %sc.mjd0 = sc.mjd0 + (tcamp/86400);                      % Updates epoch (Modified Juilan Date)
    %Ys = ode4(@(t,Y)AttDyn(t,Y,sc,sim,KOEt),[tt(i):tint:tt(i)+tcamp],Yt(i,:)');
    [~,Ys] = ode4(@(t,Y)AttDyn(t,Y,sc,sim,KOEt,jd),[tt(i):tint:tt(i)+tcamp],Yt(i,:)');
    Yt(i+1,:)=Ys(end,:);
    Yt(i+1,1:4) = Yt(i+1,1:4)/norm(Yt(i+1,1:4));
end

[~,magfield(:,1),magmoment(:,1),atterr(:,1)] = AttDyn(tt(1),Yt(1,:)',sc,sim,KOE,sc.jd0);

for j=2:length(tt)
    KOEt = kepprop2b(KOE,tt(j-1),GM);                         % Propogates orbit
    jd = sc.jd0 + ((j-1)*tcamp/86400);                        % Updates epoch (Julian Date)
    [~,magfield(:,j),magmoment(:,j),atterr(:,j)] = AttDyn(tt(j),Yt(j,:)',sc,sim,KOEt,jd);
    % Transform inertial attitude to orbital attitude
    % TBI = TBO*TOI ==> TBI*(TOI)^-1 = TBO
    % If we have inertial position and velocity from GPS, we can use a different formula
    %{
    q = qmult(Yt(j,1:4)',qinv(qmult(getq(1,-pi/2),qmult(getq(3,sc.wo*tt(i)+pi/2),getq(1,sc.inertia)))));
    att(:,i) = 2*q(1:3);
    %}
end

%% Plots
r2d = 180/pi;                                       % Converts from radians to degrees
figure
subplot(4,3,1)
plot(tt,magmoment(1,:));
title('Torque Bar 1(A*m^2)')

subplot(4,3,2)
plot(tt,magmoment(2,:));
title('Torque Bar 2(A*m^2)')

subplot(4,3,3)
plot(tt,magmoment(3,:));
title('Torque Bar 3(A*m^2)')

subplot(4,3,4)
plot([tt],Yt(:,5)');
title('Angular velocity in body-x')

subplot(4,3,5)
plot([tt],Yt(:,6)');
title('Angular velocity in body-y')

subplot(4,3,6)
plot([tt],Yt(:,7)');
title('Angular velocity in body-z')

subplot(4,3,7)
plot(tt,magfield(1,:));
title('Magnetic field component in body-x')

subplot(4,3,8)
plot(tt,magfield(2,:));
title('Magnetic field component in body-y')

subplot(4,3,9)
plot(tt,magfield(3,:));
title('Magnetic field component in body-z')

subplot(4,3,10)
plot(tt,atterr(1,:));
title('Attitude error in body-x')

subplot(4,3,11)
plot(tt,atterr(2,:));
title('Attitude error in body-y')

subplot(4,3,12)
plot(tt,atterr(3,:));
title('Attitude error in body-z')

%{
% Plots roll, pitch, yaw
figure
subplot(3,1,1)
plot(tt/60,att(1,:)*r2d)
title('Roll (deg)')
xlabel('min')

subplot(3,1,2)
plot(tt/60,att(2,:)*r2d)
title('Pitch (deg)')
xlabel('min')

subplot(3,1,3)
plot(tt/60,att(3,:)*r2d)
title('Yaw (deg)')
xlabel('min')
%}


