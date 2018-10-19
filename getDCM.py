# PRIORITY
#
"""
bV and sV converted on board
Know instantaneous time, at any measurement (possibly from the GPS)
Convert GPS time to UTC Time
Interface between data and ADCS
ecef2eci (should be from a toolbox)

"""
import numpy as np
def getDCM(bV, sV, bI, sI):
    q1 = q1.reshape(-1,1) #q1 is np.matrix
        # need to figure out how to do the convert using reshape,
        #once figured out, all of them require same steps



"""
% SHD: I recommend having a test program to test the logic of this function.
% For example, see below. I commented out so it doesn't run, but uncomment and
% paste the lines to the Command window to check. Be sure the comment the
% lines again and save this function before actually executing the line in
% the Command window.
% q = [.1 .2 -.3 1]';q=q/norm(q); % Made-up quaternion from inertial to vehicle frame
% R = q2dcm(q);
% bi = [1 2 3]';bi=bi/norm(bi); % Made-up b vector in inertial frame
% si = [-1 3 -2]';si=si/norm(si); % Made-up s vector in inertial frame
% bv = R*bi;
% sv = R*si;
% RR = getDCM(bv,sv,bi,si); % Run the made-up vectors through the routine
% dR = R-RR % This should be close to zero

"""
