function q=getq(indax,ang)
% function q=getq(indax,ang)
% Convert a rotation about a principal axis to a quaternion.
% INPUT:
%  indax  = index number of the principal axis: [1 2 3] = [x y z]
%  ang    = angle of rotation in (rad)
% OUTPUT:
%  q      = quaternion [vector;scalar]
a2=ang/2;
c=cos(a2);
s=sin(a2);
q=zeros(4,1);
q(4)=c;
if indax==1
  q(1)=s;
elseif indax==2
  q(2)=s;
elseif indax==3
  q(3)=s;
else
  error('indax must be 1, 2,or 3');
end
