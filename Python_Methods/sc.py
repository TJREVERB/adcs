from datetime import datetime

def utc2jul(dt):
    a = math.floor((14-dt.month)/12)
    y = dt.year + 4800 - a
    m = dt.month + 12*a - 3
    jdn = dt.day + math.floor((153*m + 2)/5) + 365*y + math.floor(y/4) - \
    math.floor(y/100) + math.floor(y/400) - 32045
    jd = jdn + (dt.hour - 12) / 24 + dt.minute / 1440 + \
    dt.second / 86400 - 2415020.5
    return jd

class sc:
    jd0 = utc2jul(datetime(2018, 11, 8, 12, 0, 0))
    inertia = np.diag([.0108, .0108, .0108])
