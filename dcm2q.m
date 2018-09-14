% dcm2q
% converts direction cosine matrix to quaternion [vect, scalar]

function [q] = dcm2q(dcm)

sizeDCM = size(dcm);
if(sizeDCM(1) ~= 3 || sizeDCM(2) ~= 3) 
    disp('DCM of incorrect dimensions. Quaternion initialized with dimensions of 0');
    q = [0,0,0,0];
    return;
end
q4 = .5*sqrt(1+dcm(1,1)+dcm(2,2)+dcm(3,3));
q1 = (1/(4*q4))*(dcm(2,3)-dcm(3,2));
q2 = (1/(4*q4))*(dcm(3,1)-dcm(1,3));
q3 = (1/(4*q4))*(dcm(1,2)-dcm(2,1));
q = [q1,q2,q3,q4];


