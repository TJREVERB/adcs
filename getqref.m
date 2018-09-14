% getqref
% Created by Anonto Zaman 1/26/2018
% Calculates reference quaternion (from inertial) based on keplerian 
% position such that the +Z axis points towards Earth
% Input vector kep should be formatted as follows: 
% kep = [sma ecc incl raan argp tran], where
% sma  :       Semi-major axis (m)
% ecc  :       Eccentricity
% incl :       Inclination (rad)
% raan :       Right ascension of the ascending node (rad)
% argp :       Argument of perigee (rad)
% tran :       True anomaly (rad)
%
% Updated and commented by S. Hur-Diaz 3/31/2018

% SHD: The input is really a vector, not a structure, so the 
% description above should be with respect to vector not a struct.
% Note that the angular elements are all in degrees, so they need
% to be convereted to radians before applying the trig functions
function [qref] = getqref(poskep)
% SHD: Convert angular elements from degrees to radians
% poskep(3:6) = poskep(3:6)*pi/180;
% First rotation; about the z-axis (right ascension)
q1 = [0,0,sin(poskep(4)/2),cos(poskep(4)/2)];               
% Second rotation; about the x-axis (inclination)
q2 = [sin(poskep(3)/2),0,0,cos(poskep(3)/2)];                   
% Third rotation; about the z-axis (argument of perigee and anomaly + 90 deg)
% The extra 90 deg is to align the +X axis with the velocity vector or the
% roll axis
q3 = [0,0,sin((poskep(5)+poskep(6)+(pi/2))/2),cos((poskep(5)+poskep(6)+(pi/2))/2)];   
% Fourth rotation; about the x-axis (radio facing earth)
q4 = [sin((-pi/2)/2),0,0,cos((-pi/2)/2)];      
qref = qmult(q1,q2);
qref = qmult(qref,q3);
qref = qmult(qref,q4);
qref = qref/norm(qref);