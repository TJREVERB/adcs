# Converted by Derek Goh
# October 29th, 2018

import numpy as np

<<<<<<< HEAD
"""
    #print getDCM([1,2,3,4], [1,2,3,4], [1,2,3,4], [1,2,3,4])
=======
#bV and sV converted on board
#Know instantaneous time, at any measurement (possibly from the GPS)
#Convert GPS time to UTC Time
#Interface between data and ADCS
#ecef2eci (should be from a toolbox)
>>>>>>> d3bf26f9ad2ad895c755cac24994f917a751a8b0

def normalize(v):
   norm=np.linalg.norm(v, ord=1)
   if norm==0:
       norm=np.finfo(v.dtype).eps
   return v/norm

def q2dcm(q):
   R = np.zeros((3,3))

   R[0,0] = q[0]^2-q[1]^2-q[2]^2+q[3]^2
   R[0,1] = 2*(q[0]*q[1]+q[2]*q[3])
   R[0,2] = 2*(q[0]*q[2]-q[1]*q[3])

   R[1,0] = 2*(q[0]*q[1]-q[2]*q[3])
   R[1,1] = -q[0]^2+q[1]^2-q[2]^2+q[3]^2
   R[1,2] = 2*(q[1]*q[2]+q[0]*q[3])

   R[2,0] = 2*(q[0]*q[2]+q[1]*q[3])
   R[2,1] = 2*(q[1]*q[2]-q[0]*q[3])
   R[2,2] = -q[0]^2-q[1]^2+q[2]^2+q[3]^2

   return R


def getDCM(bV, sV, bI, sI):
  bV = np.array(bV)
  sV = np.array(sV)
  bI = np.array(bI)
  sI = np.array(sI)

  bV = np.reshape(bV, -1) #
  sV = np.reshape(sV, -1) #
  bI = np.reshape(bI, -1) #
  sI = np.reshape(sI, -1) #

  vu2 = np.cross(bV, sV)
  vu2 = normalize(vu2)
  vmV = np.array(bV.getH(), vu2.getH(), np.cross(bV, vu2).getH()) #
  iu2 = np.cross(bI, sI)
  iu2 = normalize(iu2)
  imV = np.array(bI.getH(), iu2.getH(), np.cross(bI, iu2).getH()) #
  ivDCM = vmV*imV
  return ivDCM

# q = normalize([0.1 0.2 -0.3 1].getH())
# R = q2dcm(q)
# bi = normalize([1 2 3].getH())
# si = normalize([-1 3 -2].getH())
# bv = R*bi
# sv = R*si
# RR = getDCM(bv, sv, bi, si)
# dR = R - RR

print(getDCM([1,2,3],[1,2,3],[1,2,3],[1,2,3]))


# SHD: I recommend having a test program to test the logic of this function.
# For example, see below. I commented out so it doesn’t run, but uncomment and
# paste the lines to the Command window to check. Be sure to comment the
# lines again and save this function before actually executing the line in
# the Command window.
# q = [.1 .2 -.3 1]‘;q=q/norm(q); # Made-up quaternion from inertial to vehicle frame
# R = q2dcm(q);
# bi = [1 2 3]‘;bi=bi/norm(bi); # Made-up b vector in inertial frame
# si = [-1 3 -2]’;si=si/norm(si); # Made-up s vector in inertial frame
# bv = R*bi;
# sv = R*si;
# RR = getDCM(bv,sv,bi,si); # Run the made-up vectors through the routine
# dR = R-RR # This should be close to zero
