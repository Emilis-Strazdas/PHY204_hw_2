import numpy as np
from numpy import sin, cos, pi
import matplotlib.pyplot as plt
from scipy.special import roots_legendre # coefficients for Gauss-Legendre quadrature rule
from matplotlib import ticker

########## Plot formating settings ##########

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-1, 1))

########## Conditions of the system ##########

A = 1               # radius of the ring
I = 1               # current
W = 100             # angular velocity of the current
C = 10000           # speed of light

########## Integration parameters ##########

n_int = 40          # number of integration points
delta = 0.0001 / 2  # differential step size

phis, dp = roots_legendre(n_int)    # Gauss-Legendre quadrature rule
phis = (phis + 1) * pi              # Adjust phi values
dp = dp * pi                        # Adjust weights

########## Integration ##########

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
    l = A * np.array([-sin(phi), cos(phi), 0.])     # vector e in the direction of the integration
    tr = t - distance(x, y, z, phi) / C             # retarded time
    return (tr>0) * sin(W*tr) * l / distance(x, y, z, phi)

def integrate(f, x, w):
    Integral = 0

    for i in range(len(x)):
        Integral += f(x[i]) * w[i]
        
    return Integral * A * I / (C**2)

########## A - potential magnetic field ##########

def A_compute(x, y, z, t):
    return integrate(lambda phi: integrand(x, y, z, t, phi), phis, dp)

########## Partial derivatives of A ##########

def diff_A_x(x, y, z, t, d):
    return (A_compute(x + d, y, z, t) - A_compute(x - d, y, z, t)) / (2 * d)

def diff_A_y(x, y, z, t, d):
    return (A_compute(x, y + d, z, t) - A_compute(x, y - d, z, t)) / (2 * d)

def diff_A_z(x, y, z, t, d):
    return (A_compute(x, y, z + d, t) - A_compute(x, y, z - d, t)) / (2 * d)

def diff_A_t(x, y, z, t, d):
    return (A_compute(x, y, z, t + d) - A_compute(x, y, z, t - d)) / (2 * d)

########## E - electric field ##########

def E_compute(x, y, z, t, d):
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
    E = - diff_A_t(x, y, z, t, d)
    return E

########## B - magnetic field ##########

def B_compute(x, y, z, t, d):
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
    B[0] = diff_A_y(x, y, z, t, d)[2] - diff_A_z(x, y, z, t, d)[1]
    B[1] = diff_A_z(x, y, z, t, d)[0] - diff_A_x(x, y, z, t, d)[2]
    B[2] = diff_A_x(x, y, z, t, d)[1] - diff_A_y(x, y, z, t, d)[0]
    return B

########## Plot ##########

def plot(x, y, z, t0):
    t = np.linspace(t0, t0 + 8 * pi / W, 100)

    grid1 = np.zeros((len(t), 3))
    grid2 = np.zeros((len(t), 3))

    for i in range(len(t)):
        grid1[i,:] = B_compute(x, y, z, t[i], delta)
        grid2[i,:] = E_compute(x, y, z, t[i], delta)

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('$t$, $s$')

    ax1.plot(t, grid1[:, 0], label='$B_x$', color=(0/255,   0/255, 255/255), linestyle='dotted')
    ax1.plot(t, grid1[:, 1], label='$B_y$', color=(0/255, 128/255, 255/255), linestyle='dashed')
    ax1.plot(t, grid1[:, 2], label='$B_z$', color=(0/255, 255/255, 255/255), linestyle='solid')
    ax1.set_ylabel('$B$, $G$')
    ax1.yaxis.set_major_formatter(formatter)

    ax2 = ax1.twinx()
    ax2.plot(t, grid2[:, 0], label='$E_x$', color=(255/255,   0/255, 0/255), linestyle='dotted')
    ax2.plot(t, grid2[:, 1], label='$E_y$', color=(255/255, 128/255, 0/255), linestyle='dashed')
    ax2.plot(t, grid2[:, 2], label='$E_z$', color=(255/255, 255/255, 0/255), linestyle='solid')
    ax2.set_ylabel(r'$E$, $\frac{statV}{cm}$')
    ax2.yaxis.set_major_formatter(formatter)

    fig.legend()  # Adjust the legend position
    plt.title(f'Position: r = ({x}, {y}, {z}), $t_0$ = {t[0]}')

    plt.tight_layout()
    plt.show()

########## Helper functions ##########

def distance(x, y, z, phi):
    return np.sqrt( (x - A*cos(phi))**2 + (y - A*sin(phi))**2 + z**2 )

########## Main ##########

def main():
    X, Y, Z, T0 = 1000*A, 0, 0, 0
    plot(X, Y, Z, T0)

##########################

if __name__ == '__main__':
    main()