from numpy import pi
from scipy.special.orthogonal import roots_legendre
from utils.integration import *

n_int = 50          # number of integration points
delta = 0.001 / 2   # differential step size

phis, dp = roots_legendre(n_int)    # Gauss-Legendre quadrature rule
phis = (phis + 1) * pi              # Adjust phi values
dp = dp * pi                        # Adjust weights

def integrate(f):
    Int = 0

    for i in range(len(phis)):
        Int += f(phis[i]) * dp[i]
        
    return Int