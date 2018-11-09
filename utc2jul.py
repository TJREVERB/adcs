import math
from math import floor
"""
https://github.com/dannyzed/julian/blob/master/julian/julian.py
"""
def utc2jul(epoch):
    a = math.floor((14-dt.month)/12)
    y = dt.year + 4800 - a
    m = dt.month + 12*a - 3
    jdn = dt.day + math.floor((153*m + 2)/5) + 365*y + math.floor(y/4) - \
    math.floor(y/100) + math.floor(y/400) - 32045
    jd = jdn + (dt.hour - 12) / 24 + dt.minute / 1440 + \
    dt.second / 86400 - 2415020.5 #Subtracting dates since Jan 1, 1900
    # + dt.microsecond / 86400000000
    return jd
"""
a = math.floor((14-epoch[1])/12)
y = epoch[0] + 4800 - a
m = epoch[1] + 12*a - 3
jdn = epoch[2] + math.floor((153*m + 2)/5) + 365*y + math.floor(y/4) - \
math.floor(y/100) + math.floor(y/400) - 32045
jd = jdn + (epoch[3] - 12) / 24 + epoch[4] / 1440 + \
epoch[5] / 86400 - 2415020.5 #Subtracting dates since Jan 1, 1900
# + dt.microsecond / 86400000000
"""
