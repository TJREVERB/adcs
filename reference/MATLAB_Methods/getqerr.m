% getqerr
% Created by Anonto Zaman 1/26/2018
% Calculates error quaternion given a reference quaternion and a real
% quaternion
%
% Updated and commented by S. Hur-Diaz 3/31/2018

function [qerr] = getqerr(q,qref)

% Change the name of the inverse quaternion to avoid confusion
qrefinv = qinv(qref);                      % Inverse of reference quaternion

% SHD: Since q = qerr * qref, qerr = q * inv(qref). The qmult function
% expects the first argument to be the first rotation and the second
% argument to be the second rotation. The input order is opposite the
% rotation placement order in the equation. Hence, there is the confusion.
% qerr = qmult(q,qref);                   % Error quaternion is the product of the real quaternion and the inverse reference quaternion
qerr = qmult(qrefinv,q);  % Error quaternion is the product of the real quaternion and the inverse reference quaternion
qerr = qerr/norm(qerr);

% SHD: You should check the code by creating a test. For example:
% qr=[.4 .3 .1 1];qr=qr/norm(qr);
% dq=[.1 -.2 .3 1];dq=dq/norm(dq);
% q = qmult(qr,dq);
% dq_check = getqerr(q,qr);
% delta = dq-dq_check % Should be close to zero