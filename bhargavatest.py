from numpy import linalg as LA
import numpy as np
 #bV and sV converted on board
#@@ -9,12 +9,6 @@
#Interface between data and ADCS
#ecef2eci (should be from a toolbox)
def normalize(v):
   norm=np.linalg.norm(v, ord=1)
   if norm==0:
       norm=np.finfo(v.dtype).eps
   return v/norm
def q2dcm(q):
   R = np.zeros((3,3))
 #@@ -34,21 +28,21 @@ def q2dcm(q):
def getDCM(bV, sV, bI, sI):
    bV = np.array(bV)
    sV = np.array(sV)
    bI = np.array(bI)
    sI = np.array(sI)
    bV = np.matrix(bV)
    sV = np.matrix(sV)
    bI = np.matrix(bI)
    sI = np.matrix(sI)
    bV = np.reshape(bV, -1) #
    sV = np.reshape(sV, -1) #
    bI = np.reshape(bI, -1) #
    sI = np.reshape(sI, -1) #
    vu2 = np.cross(bV, sV)
    vu2 = normalize(vu2)
    vu2 = LA.norm(vu2)
    vmV = np.array(bV.getH(), vu2.getH(), np.cross(bV, vu2).getH()) #
    iu2 = np.cross(bI, sI)
    iu2 = normalize(iu2)
    iu2 = LA.norm(iu2)
    imV = np.array(bI.getH(), iu2.getH(), np.cross(bI, iu2).getH()) #
    ivDCM = vmV*imV
    return ivDCM
    #def getDCM(bV, sV, bI, sI):
# RR = getDCM(bv, sv, bi, si)
# dR = R - RR
    print(getDCM([1,2,3],[1,2,3],[1,2,3],[1,2,3]))
    #getDCM([1,2,3],[1,2,3],[1,2,3],[1,2,3])
 # SHD: I recommend having a test program to test the logic of this function.
