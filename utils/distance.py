"""
Function to calculate the distance between a point and a circle.
"""

from numpy import sin, cos, sqrt
from utils.initial_conditions import A

def distance(x, y, z, phi):
    """
    Calculates the distance between a point (x, y, z) and a point on a circle
    with radius A centered at the origin.

    Args:
        x, y, z: coordinates of the point
        phi: angle of the circle

    Returns:
        Distance between the point (x, y, z) and the point on the circle
    """
    return sqrt( (x - A * cos(phi))**2 + (y - A * sin(phi))**2 + z**2 )