% dcm2gibbs
% Created by Anonto Zaman 2/1/2018
% Conversion of direction cosine matrix to 3-1-2 gibbs vector
% Output is a row vector with components representing rotations about the
% 3-1-2 axes respectively 

function [g] = dcm2gibbs(dcm)

g = zeros([1,3]);
g(1) = -atan(dcm(2,1)/dcm(2,2));
g(2) = asin(dcm(2,3));
g(3) = -atan(dcm(1,3)/dcm(3,3));


