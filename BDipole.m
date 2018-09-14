function [b, bDot] = BDipole( r, jD, v )

%% Computes the geocentric magnetic field based on a tilted dipole model. 
% The output is in geocentric inertial coordinates (ECI). This function includes
% the effect of dipole motion on the earth.
%
% Type BDipole for a demonstration of a satellite in a low earth inclined orbit.
%--------------------------------------------------------------------------
%   Form:
%   [b, bDot] = BDipole( r, jD, v )
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   r           (3,:)   Position vector in the ECI frame (km)
%   jD          (1,:)   Julian days
%   v           (3,:)   Velocity (km/s)
%
%   -------
%   Outputs
%   -------
%   b           (3,:)   Magnetic field in the ECI frame (T)
%   bDot        (3,:)   Derivative of b in the ECI frame (T/s)
%
%--------------------------------------------------------------------------
%	Reference:  Wertz, J., ed. "Spacecraft Attitude Determination and
%               Control," Kluwer, 1976, 783-784.
%
%   Includes 1995 IGRF coefficients as of Jan. 1999
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%   Copyright (c) 1995-1999, 2014 Princeton Satellite Systems, Inc.
%   All rights reserved.
%--------------------------------------------------------------------------
%   Since version 1.
%--------------------------------------------------------------------------

% Demo
%-----
if( nargin < 2 )
    a       = 7000;
    el      = [a 55*pi/180 0 0 0 0 0];
    p       = Period(a);
    t       = linspace(0,p,1000);
    [r, v]	= RVFromKepler(el,t);
    jD      = JD2000 + t/86400;
    BDipole( r, jD, v );
    return;
end

JD2000 = juliandate([2000,1,1]);
jD1995     = -1826.5;
jDFrom2000 = jD - JD2000;

dJD       = (jDFrom2000(1)-jD1995)/365.25;

g10       = -29652   + 22.4*dJD; 
g11       =  -1787.5 + 11.3*dJD;  
h11       =   5367.5 - 15.9*dJD;  

h0        = sqrt( h11^2 + g11^2 + g10^2 ); 
a         = 6371.2;

cosThetaM = g10/h0;
thetaM    = acos( cosThetaM ); 
sinThetaM = sin( thetaM );
phiM      = atan2( h11, g11 );
uDipoleEF = [sinThetaM*cos( phiM ); sinThetaM*sin( phiM ); cosThetaM];

aCuH      = a^3*h0*1.e-9;		  
n         = length(jD);
b         = zeros(3,n);

% First path includes the derivative
%-----------------------------------
if( nargin > 2  )
    bDot      = zeros(3,n);
    for k = 1:n
        cECIToEF    = ECIToEF( JD2T( jD(k) ) );
        rEF         = cECIToEF*r(:,k);
        rMag        = Mag( rEF );
        uR          = rEF/rMag;
        bEF         = (aCuH/rMag^3)*(3*(uDipoleEF'*uR)*uR-uDipoleEF);
        b(:,k)      = cECIToEF'*bEF;
        bDot(:,k)	= DerivativeOfB( bEF, rEF, r(:,k), rMag, v(:,k), aCuH, uDipoleEF, cECIToEF', jD(k) );
    end
else
    for k = 1:n
        cECIToEF    = ECIToEF( JD2T( jD(k) ) );
        rEF         = cECIToEF*r(:,k);
        rMag        = Mag( rEF );
        uR          = rEF/rMag;
        bEF         = (aCuH/rMag^3)*(3*(uDipoleEF'*uR)*uR-uDipoleEF);
        b(:,k)      = cECIToEF'*bEF;
    end
end

% Plotting
%---------
if( nargout == 0 )
	Plot2D(jDFrom2000, r, 'Time from J2000', {'r_x (km)' 'r_y (km)' 'r_z (km)'}, 'ECI Position')
    if( nargin < 3 )
        Plot2D(jDFrom2000, b*1.e9,'Time from J2000',{'B (nT)'},'Dipole Field');
    else
        Plot2D(jDFrom2000, [b;bDot]*1.e9,'Time from J2000',{'B (nT)' 'dB/dt (nT/s)'},'Dipole Field','lin',{'[1 2 3]' '[4 5 6]'});
    end
    legend('x','y','z')
    clear b
end

%--------------------------------------------------------------------------
%   Time derivative of b
%--------------------------------------------------------------------------
function bDot = DerivativeOfB( bEF, rEF, r, rho, v, a, m, c, jD )

uEF     = rEF/rho;

rhoDot  = r'*v/rho;

omega   = [0;0;EarthRte(jD)];
cDot	  = c*Skew(omega);

uEFDot	= (c'*v + cDot'*r)/rho - rEF*rhoDot/rho^2;

bEFDot  = -3*a*rhoDot*(3*(m'*uEF)*uEF - m)/rho^4 ...
        +  3*a*((m'*uEFDot)*uEF + (m'*uEF)*uEFDot)/rho^3;
    
bDot    = c*bEFDot + cDot*bEF;


% PSS internal file version information
%--------------------------------------
% $Date: 2016-12-01 17:01:35 -0500 (Thu, 01 Dec 2016) $
% $Revision: 43640 $
