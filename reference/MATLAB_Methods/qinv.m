% qinv
% Created by Anonto Zaman 1/26/2018
% Calculates the inverse quaternion [vector,scalar]

function [q] = qinv(qin)

qin = qin(:)';                          % Converts input quaternion to row form
q = zeros([1,4]);                       % Initializes output quaternion
q(1) = -qin(1);                         % Negates vector component of input quaternion
q(2) = -qin(2);
q(3) = -qin(3);
q(4) = qin(4);

