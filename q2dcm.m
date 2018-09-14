function R=q2dcm(q)
% function R=q2dcm(q)
% Converts from a quaternion [vector;scalar] to a direction cosine matrix
R=zeros(3,3);
R(1,1)=q(1)^2-q(2)^2-q(3)^2+q(4)^2;
R(1,2)=2*(q(1)*q(2)+q(3)*q(4));
R(1,3)=2*(q(1)*q(3)-q(2)*q(4));

R(2,1)=2*(q(1)*q(2)-q(3)*q(4));
R(2,2)=-q(1)^2+q(2)^2-q(3)^2+q(4)^2;
R(2,3)=2*(q(2)*q(3)+q(1)*q(4));

R(3,1)=2*(q(1)*q(3)+q(2)*q(4));
R(3,2)=2*(q(2)*q(3)-q(1)*q(4));
R(3,3)=-q(1)^2-q(2)^2+q(3)^2+q(4)^2;
