% qmult
% Multiplies two quaternions q1 and q2 [vector, scalar]
% Created by Anonto Zaman 1/25/2018

function [q] = qmult(q1,q2)

q1 = q1(:)';
q2 = q2(:)';
comp1 = [q2(4), q2(3), -q2(2), q2(1);...
    -q2(3), q2(4), q2(1), q2(2);...
    q2(2), -q2(1), q2(4), q2(3);...
    -q2(1), -q2(2), -q2(3), q2(4)];
comp2 = [q1(1);q1(2);q1(3);q1(4)];
q = (comp1*comp2)';
