"""
Calculates the vecotr potential A at a point (x, y, z) and time t.
"""

########## Imports ##########

import numpy as np
from numpy import sin, cos
from utils.distance import distance
from utils.initial_conditions import *
from utils.integration import integrate

########## Functions ##########

def A_compute(x, y, z, t):
    """
    Calculates the vector potential A of the magnetic field B.

    Args:
        x, y, z: coordinates of the point where the vector potential is calculated
        t: time
    
    Returns:
        The vector potential A of the magnetic field B
    """
    Int = integrate(lambda phi: integrand(x, y, z, t, phi))
    return I * A / (C**2) * Int

def integrand(x, y, z, t, phi):
    """
    Integrand of the vector potential A
    
    Sets the vector e in the direction of the integration, calculates the retarded
    time and computes the integrand of the vector potential A with respect to phi.

    Args:
        x, y, z: coordinates of the point where the vector potential is calculated
        t: time
        phi: angle of the integration

    Returns:
        Integrand of the vector potential A with respect to phi, 3d vector
    """
    e = np.array([-sin(phi), cos(phi), 0.])         # vector e in the direction of the integration
    tr = t - distance(x, y, z, phi) / C             # retarded time
    return (tr>0) * sin(W*tr) * e / distance(x, y, z, phi)