# Converted by Derek Goh
# October 19th, 2018

import numpy as np
"""
bV and sV converted on board
Know instantaneous time, at any measurement (possibly from the GPS)
Convert GPS time to UTC Time
Interface between data and ADCS
ecef2eci (should be from a toolbox)

"""
<<<<<<< HEAD
def normalize(v):
    norm=np.linalg.norm(v, ord=1)
    if norm==0:
        norm=np.finfo(v.dtype).eps
    return v/norm


def getDCM(bV, sV, bI, sI):
   bV = np.matrix(bV)/normalize(bV)
   sV = np.matrix(sV)/normalize(sV)
   bI = np.matrix(bI)/normalize(bI)
   sI = np.matrix(sI)/normalize(sI)
   
   bV = np.reshape(bV, (1, -1))
   sV = np.reshape(sV, (1, -1))
   bI = np.reshape(bI, (1, -1))
   sI = np.reshape(sI, (1, -1))

   vu2 = np.cross(bV, sV)
   vu2 = vu2/normalize(vu2)
   vmV = np.matrix([bV, vu2, np.cross(bV, vu2)])
   iu2 = np.cross(bI, sI)
   iu2 = iu2/normalize(iu2)
   imV = np.matrix([bI, iu2 , np.cross(bI, iu2)])
   ivDCM = vmV*imV



# SHD: I recommend having a test program to test the logic of this function. 
# For example, see below. I commented out so it doesn't run, but uncomment and
# paste the lines to the Command window to check. Be sure the comment the
# lines again and save this function before actually executing the line in
# the Command window.
# q = [.1 .2 -.3 1]';q=q/norm(q); # Made-up quaternion from inertial to vehicle frame
# R = q2dcm(q);
# bi = [1 2 3]';bi=bi/norm(bi); # Made-up b vector in inertial frame
# si = [-1 3 -2]';si=si/norm(si); # Made-up s vector in inertial frame
# bv = R*bi;
# sv = R*si;
# RR = getDCM(bv,sv,bi,si); # Run the made-up vectors through the routine
# dR = R-RR # This should be close to zero
=======
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
>>>>>>> dce6cbf338dfb1ab5f44636076aead589d162a72
