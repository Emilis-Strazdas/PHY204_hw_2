"""
Functions to compute the derivatives of the vector potential A.
"""

########## Imports ##########
from fields.a_field import A_compute
from utils.calculus_constants import *

########## Functions ##########

def diff_A_x(x, y, z, t):
    """
    Compute the derivative of the vector potential A with respect to x.
    
    Args:
        x, y, z: coordinates of the point where the vector potential is calculated
        t: time
        
    Returns:
        float: The derivative of the vector potential A with respect to x.
    """
    return (A_compute(x + d, y, z, t) - A_compute(x - d, y, z, t)) / (2 * d)

def diff_A_y(x, y, z, t):
    """
    Compute the derivative of the vector potential A with respect to y.
    
    Args:
        x (float): x-coordinate of the point at which to compute the derivative.
        y (float): y-coordinate of the point at which to compute the derivative.
        z (float): z-coordinate of the point at which to compute the derivative.
        t (float): time at which to compute the derivative.
        
    Returns:
        float: The derivative of the vector potential A with respect to y.
    """
    return (A_compute(x, y + d, z, t) - A_compute(x, y - d, z, t)) / (2 * d)

def diff_A_z(x, y, z, t):
    """
    Compute the derivative of the vector potential A with respect to z.
    
    Args:
        x (float): x-coordinate of the point at which to compute the derivative.
        y (float): y-coordinate of the point at which to compute the derivative.
        z (float): z-coordinate of the point at which to compute the derivative.
        t (float): time at which to compute the derivative.
        
    Returns:
        float: The derivative of the vector potential A with respect to z.
    """
    return (A_compute(x, y, z + d, t) - A_compute(x, y, z - d, t)) / (2 * d)

def diff_A_t(x, y, z, t):
    """
    Compute the derivative of the vector potential A with respect to t.
    
    Args:
        x (float): x-coordinate of the point at which to compute the derivative.
        y (float): y-coordinate of the point at which to compute the derivative.
        z (float): z-coordinate of the point at which to compute the derivative.
        t (float): time at which to compute the derivative.
        
    Returns:
        float: The derivative of the vector potential A with respect to t.
    """
    return (A_compute(x, y, z, t + d) - A_compute(x, y, z, t - d)) / (2 * d)