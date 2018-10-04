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
jd = 2451545.234567
ps = jd - 2400000.5
adit = jd2gcal(2400000.5, ps)
epochvec = [[]]
for index in adit:
    epochvec[index]
print(adit[0])
