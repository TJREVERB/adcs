from datetime import datetime

def gps_get_data():
    return [{'lat': 5.3357, 'lon': -1.3483e+02, 'alt': 3.9685e+05*3.2808, 'time': datetime.utcnow(), 'x_pos': 3000000., \
        'y_pos': 5000000, 'z_pos': 1000000, 'x_vel': 3000, 'y_vel': 3000, 'z_vel': 5000}]
