import numpy as np
from numpy import linalg

def get_mc(tc,b,mmax,mtrans):
    tc = tc.getH()
    b = b.getH()                                                            
    magdip = np.cross(b,tc)/(np.linalg.norm(b)**2)
    magdip = mtrans*magdip.getH()
    magdip = magdip.getH()                                                          
    for i in range (0,len(magdip)):
        if(magdip[i]>mmax[i]):                                                      
            magdip[i] = mmax[i]
        elif(magdip[i]<(-1*mmax[i])):
            magdip[i] = -1*mmax[i]
    magdip = magdip*mtrans
    return magdip
