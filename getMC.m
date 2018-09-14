% getMC
% Created by Anonto Zaman 2/27/2018
% Calculates magnetic dipole for each torquer
% INPUTS: 
%   tc = control torque (body frame)
%   b = magnetic field (body frame)
%   mmax = maximum magnetic moment for each spacecraft magnetorquer
%   mtrans = transformation matrix from body frame to magnetic torquer frame
% OUTPUT:
%   magdip= row vector with magnetic dipole 
% Checks for magnetic saturation along each axis 

function [magdip] = getMC(tc,b,mmax,mtrans)

tc = tc(:)';                                                                      % Reformats control torque vector
b = b(:)';                                                                        % Reformats magnetic field vector

magdip = cross(b,tc)/(norm(b)^2);                                                             % Calculates magnetic dipole 
magdip = mtrans*magdip';                                                          % Transforms magnetic dipole vector into magnetic dipole axes
magdip = magdip';                                                                 % Reformats magnetic dipole into a row vector
for i=1:length(magdip)
    if(magdip(i)>mmax(i))                                                         % Checks for saturation along axes of magnetorquer
        magdip(i) = mmax(i);
    elseif(magdip(i)<-1*mmax(i))
        magdip(i) = -mmax(i); 
    end
end
magdip = magdip*mtrans;                                                      % Converts magnetic moment from magnetorquer frame to body frame