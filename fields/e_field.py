from fields.a_diffs import diff_A_t

def E_compute(x, y, z, t):
    """
    Computes the electric field E at the point (x, y, z) and time t.
    
    Using the Maxwell-Faraday equation, the electric field E is calculated as the
    negative time derivative of the vector potential A.

    Args:
        x, y, z: coordinates of the point where the electric field is calculated
        t: time
        d: differential step size
        
    Returns:
        Electric field E at the point (x, y, z) and time t, 3d vector
    """
    E = - diff_A_t(x, y, z, t)
    return E