from pyorbital import tlefile
from pyorbital.orbital import Orbital
from pyorbital.orbital import OrbitElements
from datetime import datetime

def print_header(string, key=lambda: 30):
    """
    Provides an easy way to print out a header.
    :param string: The header as a string.
    :param key: Length of the message.
    """
    print((len(string) % 2)*'-' + '{:-^{width}}'.format(string, width=key()))

def getLLA():
    tle = tlefile.read("ISS (ZARYA)", "C:\\Users\\Jason\\TJREVERB\\adcs\\orbit-predictor\\isstest.txt")
    orb = Orbital("ISS (ZARYA)")
    now = datetime.utcnow()
    print_header("Inclination:")
    print(tle.inclination)
    print_header("LLA:")
    print(orb.get_lonlatalt(now))

if __name__ == "__main__":
    getLLA()