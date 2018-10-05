"""
x = [0,1,2,3,4,5]
x = {q for q in x}
for q in x:
  print(q)
y = [q for q in x]
print(y)
QQ = {b+1:b for b in y}
print(QQ)
"""
from jdcal import gcal2jd, jd2gcal
import datetime
jd = 2458396.673843
ps = jd - 2400000.5
epochvec = list(jd2gcal(2400000.5, ps)) #Converts tuple to list
print(epochvec)

hours = int(epochvec[3]*24)
epochvec.append(epochvec[3]*24 - hours) #Sets jdtuple[4] to decimal of hours
epochvec[3] = hours

minutes = int(epochvec[4]*60)
epochvec.append(epochvec[4]*60 - minutes)
epochvec[4] = minutes

epochvec[5] = (epochvec[5]*60)
print(epochvec)
