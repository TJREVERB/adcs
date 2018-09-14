function [o1, o2, o3] = getOmega(Yt, x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%This function retrieves the angular velocity values from the simulation
o1 = Yt(5);
o2 = Yt(6); 
o3 = Yt(7);
end

