#BEING WORKED ON BY AYUSH
import numpy as np
import math
def getqref(poskep):
    for i in range(3,6):
        poskep[i]=poskep[i]*math.pi/180
    q1=np.matrix([0,0,math.sin(poskep[3]/2),math.cos(poskep[3]/2)])
    q2=np.matrix([math.sin(poskep[2]/2),0,0,math.cos(poskep[2]/2)])
    q3=np.matrix([0,0,math.sin((poskep[4]+poskep[5]+(math.pi/2))/2),math.cos((poskep[4]+poskep[5]+(math.pi/2))/2)])
    q4=np.matrix([math.sin((-math.pi/2)/2),0,0,math.cos((-math.pi/2)/2)])
    qref = qmult(q1,q2)
    qref = qmult(qref,q3)
    qref = qmult(qref,q4)
    qref = qref/normalize(qref)
    print(qref)

def normalize(v, tolerance=0.00001):
    mag2 = sum(n * n for n in v)
    if abs(mag2 - 1.0) > tolerance:
        mag = sqrt(mag2)
        v = tuple(n / mag for n in v)
    return v

def qmult(q1,q2):
    q1 = q1.reshape(-1, 1)
    q2 = q2.reshape(-1, 1)
    comp1 = np.matrix([q2[3], q2[2], -q2[1], q2[0]],[-q2[2], q2[3], q2[0], q2[1]],[q2[1], -q2[0], q2[3], q2[2]],[-q2[0], -q2[1], -q2[2], q2[3]])
    comp2 = np.matrix([q1[0]],[q1[1]],[q1[2]],[q1[3]])
    return((comp1*comp2).getH())

