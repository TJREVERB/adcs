# Written by Bharath Dileepkumar
# May 16, 2019

from . import tle_dummy

def check_tle(lat, long, alt, time):
    DEGREE_CHANGE = 5
    ALT_CHANGE = 5
    lat_in_range = lat >= -90 and lat <= 90
    lon_in_range = long >= -180 and long <= 180
    alt_in_range = alt <= 410 and alt >= 380 #Might need change
    not_in_range = not lat_in_range or not lon_in_range or not alt_in_range
    tle = tle_dummy.get_lla(time)
    if(abs(lat - tle["lat"]) > DEGREE_CHANGE or abs(long - tle["lon"]) > DEGREE_CHANGE or abs(alt - tle["alt"]) > ALT_CHANGE or not_in_range):
        return False
    else:
        return True
