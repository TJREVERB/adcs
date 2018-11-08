% getDCM
% Created by Anonto Zaman 1/19/2018
% DCM for transformation from inertial frame to vehicle frame
% Inputs are:
% bV - magnetic field vector in vehicle frame
% sV - sun vector in vehicle frame
% bI - magnetic field vector in inertial frame
% sI - sun vector in inertial frame
% All inputs are in unit vector form 
% Updated and commented by S. Hur-Diaz 3/31/2018
function [ivDCM] = getDCM(bV,sV,bI,sI)

bV = bV(:)'/norm(bV);                            % Converts all input vectors to row vectors
sV = sV(:)'/norm(sV);
bI = bI(:)'/norm(bI);
sI = sI(:)'/norm(sI);                            
vu2 = cross(bV,sV);% Cross product of magnetic field and sun vector (vehicle frame)
vu2 = vu2/norm(vu2);                    % Converts cross product into a unit vector 
vmV = [bV', vu2', cross(bV,vu2)'];      % Coordinate matrix in vehicle frame
iu2 = cross(bI,sI);                     % Cross product of magnetic field and sun vector (inertial frame)
iu2 = iu2/norm(iu2);                    % Converts cross product into a unit vector
imV = [bI', iu2', cross(bI,iu2)'];      % Coordinate matrix in inertial frame
% SHD: Fixed the formula to get the DCM from inertial frame to vehicle frame
ivDCM = vmV*imV';                        % DCM from inertial frame to vehicle frame

% SHD: I recommend having a test program to test the logic of this function. 
% For example, see below. I commented out so it doesn't run, but uncomment and
% paste the lines to the Command window to check. Be sure the comment the
% lines again and save this function before actually executing the line in
% the Command window.
% q = [.1 .2 -.3 1]';q=q/norm(q); % Made-up quaternion from inertial to vehicle frame
% R = q2dcm(q);
% bi = [1 2 3]';bi=bi/norm(bi); % Made-up b vector in inertial frame
% si = [-1 3 -2]';si=si/norm(si); % Made-up s vector in inertial frame
% bv = R*bi;
% sv = R*si;
% RR = getDCM(bv,sv,bi,si); % Run the made-up vectors through the routine
% dR = R-RR % This should be close to zero