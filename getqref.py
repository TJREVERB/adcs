#BEING WORKED ON BY AYUSH
import numpy as np
from numpy import linalg as LA
import math
def getqrefnadir(poskep):
    for i in range(2,6):
        poskep[i]=poskep[i]*math.pi/180
    q1=np.matrix([0,0,math.sin(poskep[3]/2),math.cos(poskep[3]/2)])
    q2=np.matrix([math.sin(poskep[2]/2),0,0,math.cos(poskep[2]/2)])
    q3=np.matrix([0,0,math.sin((poskep[4]+poskep[5]+(math.pi/2))/2),math.cos((poskep[4]+poskep[5]+(math.pi/2))/2)])
    q4=np.matrix([math.sin((math.pi/2)/2),0,0,math.cos((-math.pi/2)/2)])
    qref = qmult(q1,q2)
    qref = qmult(qref,q3)
    qref = qmult(qref,q4)
    qref = np.divide(qref, LA.norm(qref))
    return(qref)

def getqrefsun(poskep):
    for i in range(2,6):
        poskep[i]=poskep[i]*math.pi/180
    target=sun_vec(start_day)
    
    
    return(qref)

def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

import numpy as np
def qmult(q1,q2):
    q1 = np.squeeze(np.asarray(q1))
    q2 = np.squeeze(np.asarray(q2))
    q1 = q1.reshape(-1, 1)
    q2 = q2.reshape(-1, 1)
    comp1 = np.matrix([[q2.item(3), q2.item(2), -q2.item(1), q2.item(0)],[-q2.item(2), q2.item(3), q2.item(0), q2.item(1)],[q2.item(1), -q2.item(0), q2.item(3), q2.item(2)],[-q2.item(0), -q2.item(1), -q2.item(2), q2.item(3)]])
    comp2 = np.matrix([[q1.item(0)],[q1.item(1)],[q1.item(2)],[q1.item(3)]])
    return((comp1*comp2).getH())



