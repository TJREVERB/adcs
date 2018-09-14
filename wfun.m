% wfun
% Created by Anonto Zaman 3/25/2018
% Inputs are:
%   t = time
%   w0 = angular velocity
%   torque = total amount of torque on CubeSat
%   sc = Struct describing CubeSat
function [dwdt] = wfun(t,w0,torque,sc)

w0 = w0(:)';
torque = torque(:)';
dwdt = sc.inertia/(torque'-cross(w0',(sc.inertia*w0'))');