% getthetaerr
% Calculates error theta for each axis based on error quaternion

function [thetaerr] = getthetaerr(q)

q = q(:)';                                      % Initializes q into a row vector
thetaerr = 2*(q(1:3)/q(4));                     % Calculates attitude error in each axis by dividing q by scalar component    
