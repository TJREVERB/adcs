#!/usr/bin/env python3
import math
import sys
import numpy
import sqlite3
# sdihaspif
def juliandate(datevector):
    con = sqlite3.connect(":memory:")
    list(con.execute("select julianday(2017-01-01)"))[0][0]


def mcat(array):
    numpy.array(array)
    
def wfun(t=None, w0=None, torque=None, sc=None):
    w0 = w0(mslice[:]).cT
    torque = torque(mslice[:]).cT
    dwdt = sc.inertia / (torque.cT - cross(w0.cT, (sc.inertia * w0.cT)).cT)

def sunsensors(theta=None, phi=None, m=None):

    #Calculating the sun vector using the sun sensor with the highest voltage
    [M, I] = max(m(mslice[:]))
    # Should check if voltage is over a certain threshhold before continuing to
    # avoid calculating meaningless data when sun not in view - calibration voltage come
    # with the sun sensors

    th = theta(I)# Choose the corresponding theta value
    ph = phi(I)# Choose the corresponding phi value
    vec = mcat([sin(th) * cos(ph), OMPCSEMI, sin(th) * sin(ph), OMPCSEMI, cos(th)])#sun vector in unit form

    #For +X Panel
    if I == 1:
        R = mcat([1, 0, 0, OMPCSEMI, 0, 0, 1, OMPCSEMI, 0 - 1, 0])    #rotate the sun vector along desired axis by 90 degrees to account for sun sensor placement on cubeseat sides
        vec = R * vec

        #for -X Panel
    elif I == 2:
        R = mcat([1, 0, 0, OMPCSEMI, 0, 0 - 1, OMPCSEMI, 0, 1, 0])
        vec = R * vec

        #for +Y Panel
    elif I == 3:
        R = mcat([0, 0 - 1, OMPCSEMI, 0, 1, 0, OMPCSEMI, 1, 0, 0])
        vec = R * vec

        #for -Y Panel
    elif I == 4:
        R = mcat([0, 0, 1, OMPCSEMI, 0, 1, 0, OMPCSEMI, -1, 0, 0])
        vec = R * vec

        #No changes needed for +Z Panel
        end

        sunvec = vec    #returns sun vector in vehicle frame

def sun_vec(start_day=None):

    #  Julian days since Jan 0,1900
    #  Reference for this calculation is JD 2,415,020 which 
    #  corresponds to 12:00:00 Jan 0,1900 ET (or 12:00:00 Dec 31,1899)
    jd = 29224.5 + start_day

    #  Mean longitude of sun, measured in the ecliptic from mean 
    #  equinox of date:
    L = (279.696678 + 0.9856473354 *elmul* jd + 2.267e-13 *elmul* jd **elpow** 2)

    #  Mean anomaly of sun in radians
    Ms_r = (pi / 180) * (358.475845 + 0.985600267 *elmul* jd - (1.12e-13) *elmul* jd **elpow** 2 - (7e-20) *elmul* jd **elpow** 3)

    #  Correction between mean longitude and true longitude
    dL = 1.918 *elmul* sin(Ms_r) + 0.02 *elmul* sin(2. * Ms_r)

    #  True longitude of sun, in radians
    L_sun = rem((pi / 180) * (L + dL), 2 * pi)

    #  Compute sun unit vector in ECI frame, where the Earth's 
    #  equatorial plane is inclined inc_E radians to the ecliptic
    #  R defines a rotation about the x-axis
    inc_E = (pi / 180) * (-23.45)
    R = mcat([1, 0, 0, OMPCSEMI, 0, cos(inc_E), sin(inc_E), OMPCSEMI, 0, -sin(inc_E), cos(inc_E)])# [3,3]
    sun_ecl = mcat([cos(L_sun), OMPCSEMI, sin(L_sun), OMPCSEMI, zeros(1, size(start_day, 2))])# [3,n]
    #  Since R is constant through time, can do a simple matrix multiply
    sun_equ = R * sun_ecl# [3,n]

def qmult(q1=None, q2=None):

    q1 = q1(mslice[:]).cT
    q2 = q2(mslice[:]).cT
    comp1 = mcat([q2(4), q2(3), -q2(2), q2(1), OMPCSEMI, -q2(3), q2(4), q2(1), q2(2), OMPCSEMI, q2(2), -q2(1), q2(4), q2(3), OMPCSEMI, -q2(1), -q2(2), -q2(3), q2(4)])
    comp2 = mcat([q1(1), OMPCSEMI, q1(2), OMPCSEMI, q1(3), OMPCSEMI, q1(4)])
    q = (comp1 * comp2).cT

def qinv(qin=None):

    qin = qin(mslice[:]).cT# Converts input quaternion to row form
    q = zeros(mcat([1, 4]))# Initializes output quaternion
    q(1).lvalue = -qin(1)# Negates vector component of input quaternion
    q(2).lvalue = -qin(2)
    q(3).lvalue = -qin(3)
    q(4).lvalue = qin(4)

def q2dcm(q=None):
    # function R=q2dcm(q)
    # Converts from a quaternion [vector;scalar] to a direction cosine matrix
    R = zeros(3, 3)
    R(1, 1).lvalue = q(1) ** 2 - q(2) ** 2 - q(3) ** 2 + q(4) ** 2
    R(1, 2).lvalue = 2 * (q(1) * q(2) + q(3) * q(4))
    R(1, 3).lvalue = 2 * (q(1) * q(3) - q(2) * q(4))

    R(2, 1).lvalue = 2 * (q(1) * q(2) - q(3) * q(4))
    R(2, 2).lvalue = -q(1) ** 2 + q(2) ** 2 - q(3) ** 2 + q(4) ** 2
    R(2, 3).lvalue = 2 * (q(2) * q(3) + q(1) * q(4))

    R(3, 1).lvalue = 2 * (q(1) * q(3) + q(2) * q(4))
    R(3, 2).lvalue = 2 * (q(2) * q(3) - q(1) * q(4))
    R(3, 3).lvalue = -q(1) ** 2 - q(2) ** 2 + q(3) ** 2 + q(4) ** 2

def ode4(odefun=None, tspan=None, y0=None, *varargin):
    #ODE4  Solve differential equations with a non-adaptive method of order 4.
    #   Y = ODE4(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates 
    #   the system of differential equations y' = f(t,y) by stepping from T0 to 
    #   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
    #   The vector Y0 is the initial conditions at T0. Each row in the solution 
    #   array Y corresponds to a time specified in TSPAN.
    #
    #   Y = ODE4(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters 
    #   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...). 
    #
    #   This is a non-adaptive solver. The step sequence is determined by TSPAN
    #   but the derivative function ODEFUN is evaluated multiple times per step.
    #   The solver implements the classical Runge-Kutta method of order 4.   
    #
    #   Example 
    #         tspan = 0:0.1:20;
    #         y = ode4(@vdp1,tspan,[2 0]);  
    #         plot(tspan,y(:,1));
    #     solves the system y' = vdp1(t,y) with a constant step size of 0.1, 
    #     and plots the first component of the solution.   


    if not isnumeric(tspan):
        error(mstring('TSPAN should be a vector of integration steps.'))
        end

        if not isnumeric(y0):
            error(mstring('Y0 should be a vector of initial conditions.'))
            end

            h = diff(tspan)
            if any(sign(h(1)) * h <= 0):
                error(mstring('Entries of TSPAN are not in order.'))
                end

                # try
                f0 = feval(odefun, tspan(1), y0, varargin(mslice[:]))
                # catch
                #   msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
                #   error(msg);  
                # end  

                y0 = y0(mslice[:])            # Make a column vector.
                if not isequal(size(y0), size(f0)):
                    error(mstring('Inconsistent sizes of Y0 and f(t0,y0).'))
                    end

                    neq = length(y0)
                    N = length(tspan)
                    Y = zeros(neq, N)
                    F = zeros(neq, 4)

                    Y(mslice[:], 1).lvalue = y0
                    for i in mslice[2:N]:
                        ti = tspan(i - 1)
                        hi = h(i - 1)
                        yi = Y(mslice[:], i - 1)
                        F(mslice[:], 1).lvalue = feval(odefun, ti, yi, varargin(mslice[:]))
                        F(mslice[:], 2).lvalue = feval(odefun, ti + 0.5 * hi, yi + 0.5 * hi * F(mslice[:], 1), varargin(mslice[:]))
                        F(mslice[:], 3).lvalue = feval(odefun, ti + 0.5 * hi, yi + 0.5 * hi * F(mslice[:], 2), varargin(mslice[:]))
                        F(mslice[:], 4).lvalue = feval(odefun, tspan(i), yi + hi * F(mslice[:], 3), varargin(mslice[:]))
                        Y(mslice[:], i).lvalue = yi + (hi / 6) * (F(mslice[:], 1) + 2 * F(mslice[:], 2) + 2 * F(mslice[:], 3) + F(mslice[:], 4))
                        end
                        Y = Y.T

def kepprop2b(KOE=None, t=None, GM=None):

    # KEPPROP2B  Propagates Keplerian elements using 2-body dymanics
    #
    # INPUTS
    #   VARIABLE     SIZE   DESCRIPTION (Optional/Default)
    #      KOE       (1x1)      Structure as defined in KEPEL:
    #         .sma   (1x1)  Semi-major axis (same length units as GM)
    #         .ecc   (1x1)  Eccentricity
    #         .incl  (1x1)  Inclination (rad)
    #         .raan  (1x1)  Right ascension of the ascending node (rad)
    #         .argp  (1x1)  Argument of periapse (rad)
    #         .tran  (1x1)  True anomaly (rad)
    #      t         (1xN)      Times to propagate KOE to (secs from epoch of KOE)
    #      GM        (1x1)  Gravitational constant (same length units as .sma)
    # OUTPUTS
    #      KOEf      (1x1)  Structure as defined in KEPEL:
    #         .sma   (1xN)  Semi-major axis at times t
    #         .ecc   (1xN)  Eccentricity at times t
    #         .incl  (1xN)  Inclination at times t
    #         .raan  (1xN)  Right ascension of the ascending node at times t
    #         .argp  (1xN)  Argument of periapse at times t
    #         .tran  (1xN)  True anomaly at times t
    #
    # VALIDATION/REGRESSION TEST
    #
    #   These tests have been moved to EarthOrbitPlot_test.m to conform to
    #   the new regression testing format.
    #
    # (This file is part of ODTBX, The Orbit Determination Toolbox, and is
    #  distributed under the NASA Open Source Agreement.  See file source for
    #  more details.)

    # ODTBX: Orbit Determination Toolbox
    # 
    # Copyright (c) 2003-2011 United States Government as represented by the
    # administrator of the National Aeronautics and Space Administration. All
    # Other Rights Reserved.
    # 
    # This file is distributed "as is", without any warranty, as part of the
    # ODTBX. ODTBX is free software; you can redistribute it and/or modify it
    # under the terms of the NASA Open Source Agreement, version 1.3 or later.
    # 
    # You should have received a copy of the NASA Open Source Agreement along
    # with this program (in a file named License.txt); if not, write to the 
    # NASA Goddard Space Flight Center at opensource@gsfc.nasa.gov.

    #  REVISION HISTORY
    #   Author                  Date            Comment
    #   Kevin Berry         12/15/2008      Original
    #   Kevin Berry         06/15/2009      Corrected a bug that let the
    #                                       initial mean anomaly be ignored
    #   Kevin Berry         06/18/2009      Added self test
    #   Ravi Mathur         08/28/2012      Extracted regression test

    # These elements are being held fixed over time
    KOEf.sma = KOE.sma * ones(1, length(t))
    KOEf.ecc = KOE.ecc * ones(1, length(t))
    KOEf.incl = KOE.incl * ones(1, length(t))
    KOEf.raan = KOE.raan * ones(1, length(t))
    KOEf.argp = KOE.argp * ones(1, length(t))

    # Get the mean anomaly at the epoch of KOE
    S0.ecc = KOE.ecc
    S0.nu = KOE.tran
    M0 = kepanom(S0, mstring('M'))

    # Propagate true anomaly over times t
    T = 2 * pi * sqrt(KOE.sma ** 3 / GM)#Orbit period (sec)
    M = mod(M0 + t / T * 2 * pi, 2 * pi)#Mean anomaly (rad) at times t
    KOEf.tran = zeros(1, length(t))#initialize to increase speed

    for n in mslice[1:length(t)]:
        S.ecc = KOEf.ecc(n)
        S.M = M(n)
        KOEf.tran(n).lvalue = kepanom(S, mstring('nu'))


def kepel(r=None, v=None, GM=None):
    # KEPEL  Compute Keplerian orbital elements from Cartesian states.
    #   KOE = KEPEL(R,V) computes the two-body Keplerian orbital elements
    #   from input position vectors R and V, using the EGM-96 value of the
    #   earth's gravitational constant, GM = 3.986004415e+14 m^3/sec^2.  
    #   r and v units are meters.
    #
    #   KOE = KEPEL(R,V,GM) uses the input GM instead of the EGM-96 value.
    #   r and v units must match GM.
    #
    #   KOE = KEPEL(X,[],GM) assumes that X=[R;V].
    #   r and v units must match GM.
    #
    #   The elements are returned in the data structure KOE which has the
    #   following fields, referenced to the basis in which R & V reside:
    #      KOE.sma = semi-major axis
    #      KOE.ecc = eccentricity
    #      KOE.incl = inclination
    #      KOE.raan = right ascension of the ascending node
    #      KOE.argp = argument of periapse
    #      KOE.tran = true anomaly
    #   All angles are in radians, and SMA is in the units of R and GM.
    #
    #   If called with six output arguments instead of one, the KOE's are
    #   returned in the order listed above.
    #
    # keyword: Utilities
    #
    # See also
    #      Converting back to Cartesian states: KEP2CART
    #
    # (This file is part of ODTBX, The Orbit Determination Toolbox, and is
    #  distributed under the NASA Open Source Agreement.  See file source for
    #  more details.)

    # ODTBX: Orbit Determination Toolbox
    # 
    # Copyright (c) 2003-2011 United States Government as represented by the
    # administrator of the National Aeronautics and Space Administration. All
    # Other Rights Reserved.
    # 
    # This file is distributed "as is", without any warranty, as part of the
    # ODTBX. ODTBX is free software; you can redistribute it and/or modify it
    # under the terms of the NASA Open Source Agreement, version 1.3 or later.
    # 
    # You should have received a copy of the NASA Open Source Agreement along
    # with this program (in a file named License.txt); if not, write to the 
    # NASA Goddard Space Flight Center at opensource@gsfc.nasa.gov.

    # Steven Hughes
    # NASA GSFC

    # Russell Carpenter
    # NASA GSFC

    # Parse inputs and check for errors:
    if nargin < 3:
        GM = 3.986004415;
        print (GM)
        e + 14

        end
        [nr, nc] = size(r)
        if not xor(nr == 1, nc == 1):        # NOT(Either but not both)
            error(mstring('KEPEL:inp1NotVec'), mstring('First input must be a vector.'))
            end
            el = length(r)
            if nargin == 1 or isempty(v):
                if el != 6:
                    error(mstring('KEPEL:xNot6d'), mstring('X must be 6x1 or 1x6.'))
                    end
                    v = r(mslice[4:6])
                    r = r(mslice[1:3])
                else:
                    [nrv, ncv] = size(v)
                    if not xor(nrv == 1, ncv == 1) or length(v) != 3:
                        error(mstring('KEPEL:vNot3d'), mstring('V must be 3x1 or 1x3.'))
                        end
                        if el != 3:
                            error(mstring('KEPEL:rNot3d'), mstring('R must be 3x1 or 1x3.'))
                            end

                            # now that the lengths are checked, ensure r and v have the same
                            # orientation:
                            if nr != nrv:
                                v = v.cT                            # force v orientation to match r
                                end
                                end

                                # Compute elements (cart2oe.m code):
                                k = mcat([0, OMPCSEMI, 0, OMPCSEMI, 1])
                                h = cross(r, v)
                                n = cross(k, h)
                                N = norm(n)
                                H2 = dot(h, h)
                                V2 = dot(v, v)
                                R = norm(r)
                                e = ((V2 - GM / R) * r - (dot(r, v)) * v) / GM
                                p = H2 / GM

                                KOE.ecc = norm(e)
                                KOE.sma = p / (1 - KOE.ecc ** 2)
                                KOE.incl = acos(h(3) / sqrt(H2))
                                KOE.raan = acos(n(1) / N)
                                if n(2) < -eps:                                # Fix quadrant.
                                    KOE.raan = 2 * pi - KOE.raan
                                    end
                                    KOE.argp = acos(dot(n.cT, e) / N / KOE.ecc)
                                    if e(3) < -eps:                                    # Fix quadrant.
                                        KOE.argp = 2 * pi - KOE.argp
                                        end
                                        KOE.tran = acos(dot(e.cT, r) / KOE.ecc / R)
                                        if dot(r, v) < -eps:                                        # Fix quadrant.
                                            KOE.tran = 2 * pi - KOE.tran
                                            end
                                            KOE = orderfields(KOE, mcat([2, 1, mslice[3:6]]))
                                            if nargout == 1:
                                                varargout(1).lvalue = KOE
                                            else:
                                                varargout(1).lvalue = KOE.sma
                                                varargout(2).lvalue = KOE.ecc
                                                varargout(3).lvalue = KOE.incl
                                                varargout(4).lvalue = KOE.raan
                                                varargout(5).lvalue = KOE.argp
                                                varargout(6).lvalue = KOE.tran
                                                end

def kepanom(S=None, outstring=None):
    #  KEPANOM  Solves for one angle of anomaly given another using Kepler's Equation
    #
    #   INPUTS
    #   VARIABLE        DESCRIPTION
    #      S            Data structure containing the following fields.
    #                   Only one angle of anomaly is needed as input.
    #           S.ecc   eccentricty (unitless) -required for all conversions
    #           S.E     eccentric anomaly (radians)
    #           S.M     mean anomaly (radians)
    #           S.tran  true anomaly (radians) 
    #                   (True anomaly can also be specified in the field S.nu.)
    #
    #      outstring    String defining the desired angle of anomaly
    #                    = 'E';    % eccentric anomaly
    #                    = 'M';    % mean anomaly
    #                    = 'tran'; % true anomaly
    #                    = 'nu';   % true anomaly
    #
    #   OUTPUTS
    #      angle        Angle of the deired anomaly (eccentric, mean, or true)
    #                    in radians
    #
    #   EXAMPLES
    #      S.ecc = 0.1;
    #      S.M   = 45*pi/180;
    #      E     = kepanom(S, 'E');
    #
    # (This file is part of ODTBX, The Orbit Determination Toolbox, and is
    #  distributed under the NASA Open Source Agreement.  See file source for
    #  more details.)

    # ODTBX: Orbit Determination Toolbox
    # 
    # Copyright (c) 2003-2011 United States Government as represented by the
    # administrator of the National Aeronautics and Space Administration. All
    # Other Rights Reserved.
    # 
    # This file is distributed "as is", without any warranty, as part of the
    # ODTBX. ODTBX is free software; you can redistribute it and/or modify it
    # under the terms of the NASA Open Source Agreement, version 1.3 or later.
    # 
    # You should have received a copy of the NASA Open Source Agreement along
    # with this program (in a file named License.txt); if not, write to the 
    # NASA Goddard Space Flight Center at opensource@gsfc.nasa.gov.

    #  REVISION HISTORY
    #   Author                  Date            Comment
    #   Steven Hughes       04/09/2001      Created original subfunctions
    #   Kevin Berry         12/15/2008      Created kepanom.m based on provided
    #                                       subfunctions

    __switch_0__ = upper(outstring)
    if 0:
        pass
    elif __switch_0__ == mstring('E'):
        if any(logical_or(strcmpi(mstring('tran'), fieldnames(S)), strcmpi(mstring('nu'), fieldnames(S)))):
            angle = nu2E(S)
        elif any(strcmpi(mstring('M'), fieldnames(S))):
            angle = M2E(S)
            end

def kep2cart(*varargin):
    # KEP2CART  Convert Keplerian orbital elements to cartesian states.
    #   X = KEP2CART(KOE) converts the two-body Keplerian orbital elements in
    #   the structure KOE to Cartesian position/velocity vectors, X = [R;V],
    #   using the EGM-96 value of the earth's gravitational constant, 
    #   GM=3.986004415e+14 m^3/sec^2.  KOE is defined as in KEPEL, as follows:
    #      KOE.sma = semi-major axis
    #      KOE.ecc = eccentricity
    #      KOE.incl = inclination
    #      KOE.raan = right ascension of the ascending node
    #      KOE.argp = argument of periapse
    #      KOE.tran = true anomaly
    #   All angles are in radians, and SMA is in the length units of GM (meters
    #   by default).  The output state, X, will be expressed in the same
    #   planet-centered inertial frame to which the elements are referenced.
    #   The output will be a vector if KOE is a "scalar" structure;" otherwise
    #   it will be a matrix where each column is the position/velocity state of
    #   the corresponding set of osculating elements in the fields of KOE.
    #
    #   X = KEP2CART(KOE,GM) uses the input GM instead of the EGM-96 value.
    #
    #   Alternatively, the KOE's can be input individually in the order listed
    #   above, rather than as a structure.
    #
    # keyword: Utilities, 
    #
    # See also
    #      Converting from Cartesian states to Keplerian elements: KEPEL
    #
    # (This file is part of ODTBX, The Orbit Determination Toolbox, and is
    #  distributed under the NASA Open Source Agreement.  See file source for
    #  more details.)

    # ODTBX: Orbit Determination Toolbox
    # 
    # Copyright (c) 2003-2011 United States Government as represented by the
    # administrator of the National Aeronautics and Space Administration. All
    # Other Rights Reserved.
    # 
    # This file is distributed "as is", without any warranty, as part of the
    # ODTBX. ODTBX is free software; you can redistribute it and/or modify it
    # under the terms of the NASA Open Source Agreement, version 1.3 or later.
    # 
    # You should have received a copy of the NASA Open Source Agreement along
    # with this program (in a file named License.txt); if not, write to the 
    # NASA Goddard Space Flight Center at opensource@gsfc.nasa.gov.

    # Bob Merriam (el2xyz)
    # NASA JSC

    # Russell Carpenter
    # NASA JSC and NASA GSFC

    if nargin == 1 or nargin == 6:
        GM = 3.986004415;
        print (GM)
        e + 14

    elif nargin == 2:
        GM = varargin(2)
    else:
        GM = varargin(7)
        end
        if nargin >= 6:
            a = varargin(1)        # semi-major axis
            e = varargin(2)        # eccentricity
            i = varargin(3)        # inclination (radians)
            p = varargin(4)        # argument of perigee (radians)
            n = varargin(5)        # ascending node (radians)
            w = varargin(6)        # true anomaly (radians)
        else:
            KOE = varargin(1)






            end

            # el2xyz code (vectorized by RC):

            sp = sin(p)
            cp = cos(p)
            print (cp)

            sn = sin(n)
            cn = cos(n)
            print (cn)

            si = sin(i)
            ci = cos(i)
            print (ci)

            sw = sin(w)
            cw = cos(w)
            print (cw)


        # compute position of periapsis

            P = mcat([cn *elmul* cp - sn *elmul* sp *elmul* ci, OMPCSEMI, sn *elmul* cp + cn *elmul* sp *elmul* ci, OMPCSEMI, sp *elmul* si])

        # compute the vector at a right angle
        # to P in the direction of motion
            Q = mcat([-cn *elmul* sp - sn *elmul* cp *elmul* ci, OMPCSEMI, -sn *elmul* sp + cn *elmul* cp *elmul* ci, OMPCSEMI, cp *elmul* si])

            slr = a *elmul* (1 - e **elpow** 2)        # semi-latus rectum
            rmag = slr /eldiv/ (1 + e *elmul* cw)        # position vector magnitude

            # Position and velocity vectors:
            r(mslice[1:3], mslice[:]).lvalue = (P * diag(cw) + Q * diag(sw)) * diag(rmag)
            v(mslice[1:3], mslice[:]).lvalue = (-P * diag(sw) + Q * diag(e + cw)) * diag(sqrt(GM /eldiv/ slr))
            X = mcat([r, OMPCSEMI, v])

def getthetaerr(q=None):

    q = q(mslice[:]).cT# Initializes q into a row vector
    thetaerr = 2 * (q(mslice[1:3]) / q(4))# Calculates attitude error in each axis by dividing q by scalar component

def getqref(poskep=None):
    # SHD: Convert angular elements from degrees to radians
    # poskep(3:6) = poskep(3:6)*pi/180;
    # First rotation; about the z-axis (right ascension)
    q1 = mcat([0, 0, sin(poskep(4) / 2), cos(poskep(4) / 2)])
    # Second rotation; about the x-axis (inclination)
    q2 = mcat([sin(poskep(3) / 2), 0, 0, cos(poskep(3) / 2)])
    # Third rotation; about the z-axis (argument of perigee and anomaly + 90 deg)
    # The extra 90 deg is to align the +X axis with the velocity vector or the
    # roll axis
    q3 = mcat([0, 0, sin((poskep(5) + poskep(6) + (pi / 2)) / 2), cos((poskep(5) + poskep(6) + (pi / 2)) / 2)])
    # Fourth rotation; about the x-axis (radio facing earth)
    q4 = mcat([sin((-pi / 2) / 2), 0, 0, cos((-pi / 2) / 2)])
    qref = qmult(q1, q2)
    qref = qmult(qref, q3)
    qref = qmult(qref, q4)
    qref = qref / norm(qref)

def getqerr(q=None, qref=None):

    # Change the name of the inverse quaternion to avoid confusion
    qrefinv = qinv(qref)# Inverse of reference quaternion

    # SHD: Since q = qerr * qref, qerr = q * inv(qref). The qmult function
    # expects the first argument to be the first rotation and the second
    # argument to be the second rotation. The input order is opposite the
    # rotation placement order in the equation. Hence, there is the confusion.
    # qerr = qmult(q,qref);                   % Error quaternion is the product of the real quaternion and the inverse reference quaternion
    qerr = qmult(qrefinv, q)# Error quaternion is the product of the real quaternion and the inverse reference quaternion
    qerr = qerr / norm(qerr)

def getq(indax=None, ang=None):
    # function q=getq(indax,ang)
    # Convert a rotation about a principal axis to a quaternion.
    # INPUT:
    #  indax  = index number of the principal axis: [1 2 3] = [x y z]
    #  ang    = angle of rotation in (rad)
    # OUTPUT:
    #  q      = quaternion [vector;scalar]
    a2 = ang / 2
    c = cos(a2)
    s = sin(a2)
    q = zeros(4, 1)
    q(4).lvalue = c
    if indax == 1:
        q(1).lvalue = s
    elif indax == 2:
        q(2).lvalue = s
    elif indax == 3:
        q(3).lvalue = s
    else:
        error(mstring('indax must be 1, 2,or 3'))
        end

def getOmega(Yt=None, x=None):
    #UNTITLED Summary of this function goes here
    #   Detailed explanation goes here

    #This function retrieves the angular velocity values from the simulation
    o1 = Yt(5)
    o2 = Yt(6)
    o3 = Yt(7)
    end

def getMC1(kw=None, bdot=None, mmax=None):

    u = -kw * bdot

    for j in mslice[1:length(u)]:

        if abs(u(j)) > mmax:

            u(j).lvalue = sign(u(j)) * mmax

            end

            end

def getMC(tc=None, b=None, mmax=None, mtrans=None):

    tc = tc(mslice[:]).cT# Reformats control torque vector
    b = b(mslice[:]).cT# Reformats magnetic field vector

    magdip = cross(b, tc) / (norm(b) ** 2)# Calculates magnetic dipole
    magdip = mtrans * magdip.cT# Transforms magnetic dipole vector into magnetic dipole axes
    magdip = magdip.cT# Reformats magnetic dipole into a row vector
    for i in mslice[1:length(magdip)]:
        if (magdip(i) > mmax(i)):        # Checks for saturation along axes of magnetorquer
            magdip(i).lvalue = mmax(i)
        elif (magdip(i) < -1 * mmax(i)):
            magdip(i).lvalue = -mmax(i)
            end
            end
            magdip = magdip * mtrans        # Converts magnetic moment from magnetorquer frame to body frame

def getDCM(bV=None, sV=None, bI=None, sI=None):

    bV = bV(mslice[:]).cT / norm(bV)# Converts all input vectors to row vectors
    sV = sV(mslice[:]).cT / norm(sV)
    bI = bI(mslice[:]).cT / norm(bI)
    sI = sI(mslice[:]).cT / norm(sI)
    vu2 = cross(bV, sV)# Cross product of magnetic field and sun vector (vehicle frame)
    vu2 = vu2 / norm(vu2)# Converts cross product into a unit vector
    vmV = mcat([bV.cT, vu2.cT, cross(bV, vu2).cT])# Coordinate matrix in vehicle frame
    iu2 = cross(bI, sI)# Cross product of magnetic field and sun vector (inertial frame)
    iu2 = iu2 / norm(iu2)# Converts cross product into a unit vector
    imV = mcat([bI.cT, iu2.cT, cross(bI, iu2).cT])# Coordinate matrix in inertial frame
    # SHD: Fixed the formula to get the DCM from inertial frame to vehicle frame
    ivDCM = vmV * imV.cT# DCM from inertial frame to vehicle frame


bv = mcat([0, 0, 0])# Initializes empty vector for placeholder magnetic field

def dcm2q(dcm=None):

    sizeDCM = size(dcm)
    if (sizeDCM(1) != 3 or sizeDCM(2) != 3):
        disp(mstring('DCM of incorrect dimensions. Quaternion initialized with dimensions of 0'))
        q = mcat([0, 0, 0, 0])
        return
        end
        q4 = .5 * sqrt(1 + dcm(1, 1) + dcm(2, 2) + dcm(3, 3))
        q1 = (1 / (4 * q4)) * (dcm(2, 3) - dcm(3, 2))
        q2 = (1 / (4 * q4)) * (dcm(3, 1) - dcm(1, 3))
        q3 = (1 / (4 * q4)) * (dcm(1, 2) - dcm(2, 1))
        q = mcat([q1, q2, q3, q4])

def dcm2gibbs(dcm=None):

    g = zeros(mcat([1, 3]))
    g(1).lvalue = -atan(dcm(2, 1) / dcm(2, 2))
    g(2).lvalue = asin(dcm(2, 3))
    g(3).lvalue = -atan(dcm(1, 3) / dcm(3, 3))

def crossm(u=None):
    # function M=crossm(u)
    # Computes the cross product matrix of the vector u


    M = mcat([0 - u(3), u(2), OMPCSEMI, u(3), 0 - u(1), OMPCSEMI, -u(2), u(1), 0])

def BDipole(r=None, jD=None, v=None):

    #% Computes the geocentric magnetic field based on a tilted dipole model. 
    # The output is in geocentric inertial coordinates (ECI). This function includes
    # the effect of dipole motion on the earth.
    #
    # Type BDipole for a demonstration of a satellite in a low earth inclined orbit.
    #--------------------------------------------------------------------------
    #   Form:
    #   [b, bDot] = BDipole( r, jD, v )
    #--------------------------------------------------------------------------
    #
    #   ------
    #   Inputs
    #   ------
    #   r           (3,:)   Position vector in the ECI frame (km)
    #   jD          (1,:)   Julian days
    #   v           (3,:)   Velocity (km/s)
    #
    #   -------
    #   Outputs
    #   -------
    #   b           (3,:)   Magnetic field in the ECI frame (T)
    #   bDot        (3,:)   Derivative of b in the ECI frame (T/s)
    #
    #--------------------------------------------------------------------------
    #   Reference:  Wertz, J., ed. "Spacecraft Attitude Determination and
    #               Control," Kluwer, 1976, 783-784.
    #
    #   Includes 1995 IGRF coefficients as of Jan. 1999
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    #   Copyright (c) 1995-1999, 2014 Princeton Satellite Systems, Inc.
    #   All rights reserved.
    #--------------------------------------------------------------------------
    #   Since version 1.
    #--------------------------------------------------------------------------

    # Demo
    #-----
    if (nargin < 2):
        a = 7000
        el = mcat([a, 55 * pi / 180, 0, 0, 0, 0, 0])
        p = Period(a)
        t = linspace(0, p, 1000)
        [r, v] = RVFromKepler(el, t)
        jD = JD2000 + t / 86400
        BDipole(r, jD, v)
        return
    end

    JD2000 = juliandate(mcat([2000, 1, 1]))
    jD1995 = -1826.5
    jDFrom2000 = jD - JD2000

    dJD = (jDFrom2000(1) - jD1995) / 365.25

    g10 = -29652 + 22.4 * dJD
    g11 = -1787.5 + 11.3 * dJD
    h11 = 5367.5 - 15.9 * dJD

    h0 = sqrt(h11 ** 2 + g11 ** 2 + g10 ** 2)
    a = 6371.2

    cosThetaM = g10 / h0
    thetaM = acos(cosThetaM)
    sinThetaM = sin(thetaM)
    phiM = atan2(h11, g11)
    uDipoleEF = mcat([sinThetaM * cos(phiM), OMPCSEMI, sinThetaM * sin(phiM), OMPCSEMI, cosThetaM])

    aCuH = a ** 3 * h0 * 1.e-9
    n = length(jD)
    b = zeros(3, n)

    # First path includes the derivative
    #-----------------------------------
    if (nargin > 2):
        bDot = zeros(3, n)
        for k in mslice[1:n]:
            cECIToEF = ECIToEF(JD2T(jD(k)))
            rEF = cECIToEF * r(mslice[:], k)
            rMag = Mag(rEF)
            uR = rEF / rMag
            bEF = (aCuH / rMag ** 3) * (3 * (uDipoleEF.cT * uR) * uR - uDipoleEF)
            b(mslice[:], k).lvalue = cECIToEF.cT * bEF
            bDot(mslice[:], k).lvalue = DerivativeOfB(bEF, rEF, r(mslice[:], k), rMag, v(mslice[:], k), aCuH, uDipoleEF, cECIToEF.cT, jD(k))
        end
    else:
        for k in mslice[1:n]:
            cECIToEF = ECIToEF(JD2T(jD(k)))
            rEF = cECIToEF * r(mslice[:], k)
            rMag = Mag(rEF)
            uR = rEF / rMag
            bEF = (aCuH / rMag ** 3) * (3 * (uDipoleEF.cT * uR) * uR - uDipoleEF)
            b(mslice[:], k).lvalue = cECIToEF.cT * bEF
        end
    end

    # Plotting
    #---------
    if (nargout == 0):
        Plot2D(jDFrom2000, r, mstring('Time from J2000'), mcellarray([mstring('r_x (km)'), mstring('r_y (km)'), mstring('r_z (km)')]), mstring('ECI Position'))
        if (nargin < 3):
            Plot2D(jDFrom2000, b * 1.e9, mstring('Time from J2000'), mcellarray([mstring('B (nT)')]), mstring('Dipole Field'))
        else:
            Plot2D(jDFrom2000, mcat([b, OMPCSEMI, bDot]) * 1.e9, mstring('Time from J2000'), mcellarray([mstring('B (nT)'), mstring('dB/dt (nT/s)')]), mstring('Dipole Field'), mstring('lin'), mcellarray([mstring('[1 2 3]'), mstring('[4 5 6]')]))
        end
        legend(mstring('x'), mstring('y'), mstring('z'))
        clear(mstring('b'))
    end

    #--------------------------------------------------------------------------
    #   Time derivative of b
    #--------------------------------------------------------------------------
#@mfunction("bDot")
def DerivativeOfB(bEF=None, rEF=None, r=None, rho=None, v=None, a=None, m=None, c=None, jD=None):

    uEF = rEF / rho

    rhoDot = r.cT * v / rho

    omega = mcat([0, OMPCSEMI, 0, OMPCSEMI, EarthRte(jD)])
    cDot = c * Skew(omega)

    uEFDot = (c.cT * v + cDot.cT * r) / rho - rEF * rhoDot / rho ** 2

    bEFDot = -3 * a * rhoDot * (3 * (m.cT * uEF) * uEF - m) / rho ** 4 + 3 * a * ((m.cT * uEFDot) * uEF + (m.cT * uEF) * uEFDot) / rho ** 3

    bDot = c * bEFDot + cDot * bEF


    # PSS internal file version information
    #--------------------------------------
    # $Date: 2016-12-01 17:01:35 -0500 (Thu, 01 Dec 2016) $
    # $Revision: 43640 $

# AttDyn
# Created by Anonto Zaman  3/25/2018
# Function that describes the differential equation of q
# Inputs are: 
#   t = time
#   Y = input vector in the form [q, w] where q represents quaternion in
#   vehicle frame and w represents angular velocity in vehicle frame
#   sc = struct of spacecraft properties
#   sim = struct of simulation properties(gain, maximum magnetic moment,
#   transformation matrix from body to magnetorquer frame)
#   KOE = struct of spacecraft keplerian orbital elements
#@mfunction("dydt, bV, magdip, thetaerr")
def AttDyn(t=None, Y=None, sc=None, sim=None, KOE=None, jd=None):

    Y = Y(mslice[:]).cT# Converts Y to row vector
    q0 = Y(mslice[1:4])# Quaternion rotation from inertial to vehicle frame
    w0 = Y(mslice[5:7])# Angular velocity vector of satellite in vehicle frame
    dcm = q2dcm(q0)# DCM from inertial to vehicle frame

    # magdip
    # SHD: Move epoch to main calling function

    #% Magnetic Field Model
    epochvec = datevec(datetime(jd, mstring('convertfrom'), mstring('juliandate')))
    cart = kep2cart(KOE)# Cartesian coordinates of satellite [pos,vel]
    cartloc = cart(mslice[1:3])# Cartesian location of satellite (row vector form)
    lla = eci2lla(cartloc.cT, epochvec)
    eci2ecef = dcmeci2ecef(mstring('IAU-2000/2006'), epochvec)
    ecef2eci = eci2ecef.cT
    magECEF = wrldmagm(lla(3), lla(1), lla(2), decyear(mstring('01-January-2018'), mstring('dd-mm-yyyy')))
    magECI = ecef2eci * magECEF
    bV = 1.0e-09 * dcm * magECI

    #% Control Torque Calculation
    #qref = getqref(struct2array(KOE));
    qref = mcat([0, 0, -sqrt(2) / 2, sqrt(2) / 2])
    qerr = getqerr(q0, qref)
    thetaerr = getthetaerr(qerr)
    #ctcomm = -1*sim.gain*thetaerr-1*sim.rgain*(w0-sim.wref);                                  % Expected control torque (row vector)
    ctcomm = -1 * sim.gain * thetaerr
    magdip = getMC(ctcomm.cT, bV, sim.mmax, sim.mtrans)# Magnetic dipole in body frame
    ctprod = cross(magdip, bV)

    #% dydt component
    dydt(mslice[1:4]).lvalue = .5 * qmult(mcat([w0, 0]), q0)# Initializes first part of q
    dydt(mslice[5:7]).lvalue = sc.inertia(ctprod.cT - cross(w0.cT, (sc.inertia * w0.cT)))
    dydt = dydt.cT

# Assumes no disturbance torques, no gravity gradient 

#% Spacecraft properties
epoch = mcat([2018, 1, 1, 12, 00, 00])# Time in UTC (3/17/2018 12:00:00 UTC)mjd = mjuliandate([y,mo,d,h,mi,s])
sc.jd0 = juliandate(epoch)# Initial satellite Julian Date
#sc.mjd0 = mjuliandate(epoch);                        % Iniitial Modified Julian Date
sc.inertia = diag(mcat([.0108, .0108, .00426]))# Assuming mass of 2.6kg and side length of .1x.1x.2 meters
# sc.rc = 6771;                                        % Orbit radius (km) (LEO), based on ISS average altitude
#sc.In = 51.6*pi/180;                                 % Orbit inclination (rad), based on ISS inclination
#sc.Torb = 2*pi/sqrt(3.986e5/sc.rc^3);                % Orbital period (s), based on calculated orbital period from ISS altitude, note decay

#% Orbital properties
GM = 3.986004418e14# (m^3/s^2)
KOE.sma = 6778557# Semi-major axis (m)
KOE.ecc = 0.0001973# Eccentricity
KOE.incl = 51.6397 * pi / 180# Inclination (rad)
KOE.raan = 115.7654 * pi / 180# Right ascension of the ascending node (rad)
KOE.argp = 211.4532 * pi / 180# Argument of perigee (rad)
KOE.tran = 0 * pi / 180# True anomaly (rad)
koeVect = struct2array(KOE)# Keplerian orbital elements in vector form
epoch = mcat([2018, 1, 1, 12, 00, 00])# Time in UTC (3/17/2018 12:00:00 UTC)
cart = kep2cart(KOE)# Cartesian coordinates of satellite [pos,vel]
cartloc = cart(mslice[1:3])# Cartesian location of satellite (row vector form)
#lla = eci2lla(cartloc',epoch);                       % Converts cartesian (ECI) coordinates to Latitude, Longitude, Altitude coordinates

#% Magnetic Field Model
epochvec = datevec(datetime(sc.jd0, mstring('convertfrom'), mstring('juliandate')))
lla = eci2lla(cartloc.cT, epochvec)
eci2ecef = dcmeci2ecef(mstring('IAU-2000/2006'), epochvec)
ecef2eci = eci2ecef.cT
magECEF = wrldmagm(lla(3), lla(1), lla(2), decyear(mstring('01-January-2018'), mstring('dd-mm-yyyy')))
magECI = ecef2eci * magECEF

#% Initial CubeSat Attitude
#qtrue = [.5,.5,.5,.99];                     % Initial rotation quaternion (rotation from Inertial to Vehicle Frame)
#qtrue = [0.5435   -0.0028   -0.6124   -0.5741];
qtrue = mcat([0, 0, sqrt(2) / 2, sqrt(2) / 2])
qtrue = qtrue / norm(qtrue)# Formats rotation quaternion
DCMtrue = q2dcm(qtrue)# Initial true quaternion in DCM form

#% Sensor Outputs
#[magTotal,~] = BDipole(cart,sc.jd0,[0;0;0]);        % Magnetic field in nT (inputs are location in ECI frame, Julian Date, and Initial Angular rate)
bI = 1.0e-09 * magECI# Magnetic field vector in ECI frame (measured in T)
bI = bI / norm(bI)
sI = sun_vec(sc.jd0 - juliandate(mcat([1980, 1, 6])))# Sun vector in ECI frame
sI = sI / norm(sI)
bV = DCMtrue * bI# Magnetic field vector in vehicle frame
bV = bV / norm(bV)
sV = DCMtrue * sI# Sun vector in vehicle frame
sV = sV / norm(sV)

#% Attitude properties
dcm = getDCM(bV, sV, bI, sI)# DCM from vehicle to inertial frame
q = dcm2q(dcm)# Spacecraft initial quaternion
#{
qref = getqref(koeVect)# Calculates reference quaternion
qerr = getqerr(q, qref)# Calcualtes error quaternion
thetaerr = getthetaerr(qerr)# Error angles for each axis
#}

#% Simulation parameters
Ts = 100# Simulation time (sec)
tcamp = 1# Plot sampling time/integration step
tint = .1# Integration step
tt = mslice[0:tcamp:Ts]# Time vector for simulation

#% Initial angular rates
n = sqrt(GM / (KOE.sma ** 3))
#w = [0,-n,0];
w = mcat([0, 0, 0])
Yt = mcat([q, w, OMPCSEMI, zeros(length(tt) - 1, 7)])

#% Control law 
sim.gain = 2e-5# Test gain value
sim.rgain = 2e-5
sim.mmax = mcat([.2, .2, .2])# Maximum magnetic moment for each torquer axis
sim.mtrans = mcat([1, 0, 0, OMPCSEMI, 0, 1, 0, OMPCSEMI, 0, 0, 1])# DCM to convert from body frame to magnetorquer frame
sim.wref = w
#{
magdip = getMC(ctcomm, bV, mmax, mtrans)# Magnetic dipole in body frame
ctprod = cross(magdip, bV)# Produced control torque
#}

#% Inputs for plotting
magmoment = zeros(3, length(tt))
magfield = zeros(3, length(tt))
atterr = zeros(3, length(tt))

#% Simulation

for i in mslice[1:length(tt) - 1]:
    KOEt = kepprop2b(KOE, tt(i), GM)# Propogates orbit
    jd = sc.jd0 + (i * tcamp / 86400)# Updates epoch (Julian Date)
    #sc.mjd0 = sc.mjd0 + (tcamp/86400);                      % Updates epoch (Modified Juilan Date)
    #Ys = ode4(@(t,Y)AttDyn(t,Y,sc,sim,KOEt),[tt(i):tint:tt(i)+tcamp],Yt(i,:)');
    ode4(lambda t, Y: AttDyn(t, Y, sc, sim, KOEt, jd), mcat([mslice[tt(i):tint:tt(i) + tcamp]]), Yt(i, mslice[:]).cT)
    Yt(i + 1, mslice[:]).lvalue = Ys(end, mslice[:])
    Yt(i + 1, mslice[1:4]).lvalue = Yt(i + 1, mslice[1:4]) / norm(Yt(i + 1, mslice[1:4]))
end

magfield(mslice[:], 1)
magmoment(mslice[:], 1)
AttDyn(tt(1), Yt(1, mslice[:]).cT, sc, sim, KOE, sc.jd0)

for j in mslice[2:length(tt)]:
    KOEt = kepprop2b(KOE, tt(j - 1), GM)# Propogates orbit
    jd = sc.jd0 + ((j - 1) * tcamp / 86400)# Updates epoch (Julian Date)
    magfield(mslice[:], j)
    magmoment(mslice[:], j)
    AttDyn(tt(j), Yt(j, mslice[:]).cT, sc, sim, KOEt, jd)
    # Transform inertial attitude to orbital attitude
    # TBI = TBO*TOI ==> TBI*(TOI)^-1 = TBO
    # If we have inertial position and velocity from GPS, we can use a different formula
    #{
    q = qmult(Yt(j, mslice[1:4]).cT, qinv(qmult(getq(1, -pi / 2), qmult(getq(3, sc.wo * tt(i) + pi / 2), getq(1, sc.inertia)))))
    att(mslice[:], i).lvalue = 2 * q(mslice[1:3])
    #}
end

#% Plots
r2d = 180 / pi# Converts from radians to degrees
figure()
subplot(4, 3, 1)
plot(tt, magmoment(1, mslice[:]))
title(mstring('Torque Bar 1(A*m^2)'))

subplot(4, 3, 2)
plot(tt, magmoment(2, mslice[:]))
title(mstring('Torque Bar 2(A*m^2)'))

subplot(4, 3, 3)
plot(tt, magmoment(3, mslice[:]))
title(mstring('Torque Bar 3(A*m^2)'))

subplot(4, 3, 4)
plot(mcat([tt]), Yt(mslice[:], 5).cT)
title(mstring('Angular velocity in body-x'))

subplot(4, 3, 5)
plot(mcat([tt]), Yt(mslice[:], 6).cT)
title(mstring('Angular velocity in body-y'))

subplot(4, 3, 6)
plot(mcat([tt]), Yt(mslice[:], 7).cT)
title(mstring('Angular velocity in body-z'))

subplot(4, 3, 7)
plot(tt, magfield(1, mslice[:]))
title(mstring('Magnetic field component in body-x'))

subplot(4, 3, 8)
plot(tt, magfield(2, mslice[:]))
title(mstring('Magnetic field component in body-y'))

subplot(4, 3, 9)
plot(tt, magfield(3, mslice[:]))
title(mstring('Magnetic field component in body-z'))

subplot(4, 3, 10)
plot(tt, atterr(1, mslice[:]))
title(mstring('Attitude error in body-x'))

subplot(4, 3, 11)
plot(tt, atterr(2, mslice[:]))
title(mstring('Attitude error in body-y'))

subplot(4, 3, 12)
plot(tt, atterr(3, mslice[:]))
title(mstring('Attitude error in body-z'))

#{
# Plots roll, pitch, yaw
figure()
subplot(3, 1, 1)
plot(tt / 60, att(1, mslice[:]) * r2d)
title(mstring('Roll (deg)'))
xlabel(mstring('min'))

subplot(3, 1, 2)
plot(tt / 60, att(2, mslice[:]) * r2d)
title(mstring('Pitch (deg)'))
xlabel(mstring('min'))

subplot(3, 1, 3)
plot(tt / 60, att(3, mslice[:]) * r2d)
title(mstring('Yaw (deg)'))
xlabel(mstring('min'))
#}
