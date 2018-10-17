#BEING WORKED ON BY AYUSH
import numpy as np
def qmult(q1,q2):
    q1 = q1.reshape(-1, 1)
    q2 = q2.reshape(-1, 1)
    comp1 = np.matrix([q2(4), q2(3), -q2(2), q2(1)],[-q2(3), q2(4), q2(1), q2(2)],[q2(2), -q2(1), q2(4), q2(3)],[-q2(1), -q2(2), -q2(3), q2(4)])
    comp2 = np.matrix([q1(1)],[q1(2)],[q1(3)],[q1(4)])
    return((comp1*comp2).getH())

