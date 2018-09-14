% qfun
% Created by Anonto Zaman  3/25/2018
% Function that describes the differential equation of q
% Inputs are: 
%   t = time
%   q0 = initial quaternion in vehicle frame
%   w0 = angular velocity in vehicle frame
function [dydt] = AttDyn(t,Y,w0)

w0 = w0(:)';
dqdt = .5*qmult([w0,0],q0);
