from numpy import sin, cos, sqrt

def distance(x, y, z, phi):
    return sqrt( (x - cos(phi))**2 + (y - sin(phi))**2 + z**2 )