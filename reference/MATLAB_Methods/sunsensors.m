% sunsensors
%
% Given angle of incidences for all 5 sun sensors on 2U cubesat
%
% The product information sheet for the ISIS sun sensors we are using 
% states that it outputs analog channels from the 4 quadrant
% photovoltic reciever, but also states that the sensor comes with an algorithm to
% compute the angle of incidence from those voltages and can give us calibration
% voltages. Since we don't yet have access to the algorithm and the
% expected voltage ranges and calibration data that come with the sensors,
% I've created a function to determine the 2 angles of incidence
% given the 4 voltages and return these angles of incidence in
% inp2ang.m, but the actual values will differ depending on the calibration voltages
% and specifics tht come with the sun sensors.
%
% Inputs:
%   theta = input vector with angle of incidence on the ZY plane for each sun sensor
%   phi = input vector with angle of incidence on the XY plane for each sun sensor
%   m = input vector with magnitude (voltage on 4 quadrants) for each sun vector
%   input vector order [+X panel, -X panel, +Y panel, -Y panel, +Z panel]

function[sunvec] = sunsensors(theta, phi, m)

%Calculating the sun vector using the sun sensor with the highest voltage
[M,I]=max(m(:)); 
% Should check if voltage is over a certain threshhold before continuing to
% avoid calculating meaningless data when sun not in view - calibration voltage come
% with the sun sensors

th = theta(I); % Choose the corresponding theta value
ph = phi(I); % Choose the corresponding phi value
vec = [sin(th)*cos(ph);sin(th)*sin(ph);cos(th)]; %sun vector in unit form

%For +X Panel
if I == 1
    R = [1 0 0; 0 0 1; 0 -1 0]; %rotate the sun vector along desired axis by 90 degrees to account for sun sensor placement on cubeseat sides
    vec = R*vec;
    
%for -X Panel
elseif I == 2
    R = [1 0 0; 0 0 -1; 0 1 0];
    vec = R*vec;
    
%for +Y Panel
elseif I == 3
    R = [0 0 -1; 0 1 0; 1 0 0];
    vec = R*vec;
    
%for -Y Panel
elseif I == 4
    R = [0 0 1; 0 1 0; -1 0 0];
    vec = R*vec;
    
%No changes needed for +Z Panel
end

sunvec = vec; %returns sun vector in vehicle frame


