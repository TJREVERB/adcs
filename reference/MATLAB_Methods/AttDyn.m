% AttDyn
% Created by Anonto Zaman  3/25/2018
% Function that describes the differential equation of q
% Inputs are: 
%   t = time
%   Y = input vector in the form [q, w] where q represents quaternion in
%   vehicle frame and w represents angular velocity in vehicle frame
%   sc = struct of spacecraft properties
%   sim = struct of simulation properties(gain, maximum magnetic moment,
%   transformation matrix from body to magnetorquer frame)
%   KOE = struct of spacecraft keplerian orbital elements
%   
function [dydt,bV,magdip,thetaerr] = AttDyn(t,Y,sc,sim,KOE,jd)

Y = Y(:)';                                      % Converts Y to row vector
q0 = Y(1:4);                                    % Quaternion rotation from inertial to vehicle frame
w0 = Y(5:7);                                    % Angular velocity vector of satellite in vehicle frame
dcm = q2dcm(q0);                                % DCM from inertial to vehicle frame 

% magdip
% SHD: Move epoch to main calling function

%% Magnetic Field Model
epochvec = datevec(datetime(jd,'convertfrom','juliandate'));
cart = kep2cart(KOE);                                % Cartesian coordinates of satellite [pos,vel]
cartloc = cart(1:3);                                 % Cartesian location of satellite (row vector form)
lla = eci2lla(cartloc',epochvec); 
eci2ecef = dcmeci2ecef('IAU-2000/2006',epochvec);
ecef2eci = eci2ecef';
magECEF = wrldmagm(lla(3),lla(1),lla(2),decyear('01-January-2018','dd-mm-yyyy'));
magECI = ecef2eci*magECEF;  
bV = 1.0e-09*dcm* magECI;

%% Control Torque Calculation
%qref = getqref(struct2array(KOE));
qref = [0,0,-sqrt(2)/2,sqrt(2)/2];
qerr = getqerr(q0,qref);
thetaerr = getthetaerr(qerr);
%ctcomm = -1*sim.gain*thetaerr-1*sim.rgain*(w0-sim.wref);                                  % Expected control torque (row vector)
ctcomm = -1*sim.gain*thetaerr;
magdip = getMC(ctcomm',bV,sim.mmax,sim.mtrans);               % Magnetic dipole in body frame
ctprod = cross(magdip,bV);

%% dydt component
dydt(1:4) = .5*qmult([w0,0],q0);                % Initializes first part of q  
dydt(5:7) = sc.inertia\(ctprod'-cross(w0',(sc.inertia*w0')));
dydt = dydt';

