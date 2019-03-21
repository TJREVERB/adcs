function M=crossm(u)
% function M=crossm(u)
% Computes the cross product matrix of the vector u
M=[0 -u(3) u(2);
   u(3) 0 -u(1);
   -u(2) u(1) 0];

    