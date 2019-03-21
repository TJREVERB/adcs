% INPUTS: cartesian state vectors:   
% Position: r(t)
% Velocity: dr/dt
% Gravitational constant
% Mass of earth
% 
% OUTPUTS: Keplarian elements [a escalar w l i m]
% Semi-major axis a
% Eccentricity e
% Argument of periapsis(perigee) w
% Longitude of ascending node l
% Inclination i
% Mean anomaly m
% Epoch is not detemined in this function as of yet
function [kepl]=cart2kep(r, vel)
    G = 6.67408 * 10^(-11); %Gravitational constant(N kg^-2 m^2)
    M = 5.9722 * 10^(24); %Mass of earth(kg)
    GM = G*M; %Constant parameter
    k = [0 0 1]; %Constant(Z-axis) used with angular momentum to calculate ascending node
    h = cross(r,vel);%Orbital momentum(Without mass dimenstion): r cross vel
    n = cross(k, h); %vector from Earth center to ascending node
    
    %Eccentricity vector and eccentricity: (e, escalar)
    e = cross(vel, h)/GM - r/norm(r);
    escalar = norm(e);
    
    %True anomaly: (v)
    if dot(r, vel) >= 0
        v = acos(dot(e,r)/(norm(e) * norm(r)));
    else
        v = 2*pi - acos(dot(e,r)/(norm(e) * norm(r)));
    end
    
    %Inclination: (i)
    i = acos(h(3)/norm(h));
    
    %Eccentric anomaly: (E) if u don't understand: https://en.wikipedia.org/wiki/Eccentric_anomaly
    E = 2*atan(tan(v/2)/sqrt((1+escalar)/(1-escalar)));
    
    %Longitude of the ascending node: (l)
    if n(2) >= 0
        l = acos(n(1)/norm(n));
    else
        l = 2*pi - acos(n(1)/norm(n));
    end
    
    %Argument of periapsis(Perigee) (w)
    if e(3) >= 0
        w = acos(dot(n,e)/(norm(n) * norm(e)));
    else
        w = 2*pi - acos(dot(n,e)/(norm(n) * norm(e)));
    end
    
    %Mean anomaly: (m)
    m = E - escalar*sin(E);
    
    %Semi-major axis: (a)
    a = 1/(2/norm(r) - (norm(vel)^2)/GM);
    
    kepl = [a escalar w l i m];