"""
Calculates the magnetic field B at a point (x, y, z) and time t.
"""

########## Imports ##########

import numpy as np
from fields.a_diffs import *

########## Functions ##########

def B_compute(x, y, z, t):
    """
    Computes the magnetic field B at the point (x, y, z) and time t.
    
    Using the Maxwell-Ampere equation, the magnetic field B is calculated as the
    curl of the vector potential A.
    
    Args:
        x, y, z: coordinates of the point where the magnetic field is calculated
        t: time
        d: differential step size
        
    Returns:
        Magnetic field B at the point (x, y, z) and time t, 3d vector
    """
    B = np.zeros(3)
    B[0] = diff_A_y(x, y, z, t)[2] - diff_A_z(x, y, z, t)[1]
    B[1] = diff_A_z(x, y, z, t)[0] - diff_A_x(x, y, z, t)[2]
    B[2] = diff_A_x(x, y, z, t)[1] - diff_A_y(x, y, z, t)[0]
    return B