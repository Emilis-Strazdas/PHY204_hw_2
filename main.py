import numpy as np
from numpy import sin, cos, pi
import matplotlib.pyplot as plt
from scipy.special import roots_legendre # coefficients for Gauss-Legendre quadrature rule
from init_conditions import *

########## Conditions of the system ##########

a = 1               # radius of the ring
I = 1               # current
w = 100             # angular velocity of the current
c = 10000           # speed of light

########## Integration parameters ##########

n_int = 50          # number of integration points
delta = 0.001 / 2   # differential step size

phis, dp = roots_legendre(n_int)    # Gauss-Legendre quadrature rule
phis = (phis + 1) * pi / 2          # phi values

########## Integration ##########

def intergrand(x, y, z, t, phi):
    e = np.array([-sin(phi), cos(phi), 0.])
    tr = t - distance(x, y, z, phi) / c

    if tr < 0:
        return np.zeros(3)
    
    return (tr>0) * sin(w*tr) * e * I * a / distance(x, y, z, phi)

def integrate(f, x, w):
    Integral = 0

    for i in range(len(x)):
        Integral += f(x[i]) * w[i]

    return Integral

########## A - potential magnetic field ##########

def A_compute(x, y, z, t):
    return integrate(lambda phi: intergrand(x, y, z, t, phi), phis, dp) / c**2
    return integrate(lambda phi: intergrand(x, y, z, t, phi), phis, dp)

def Afield(x, y, z, t):
    nphi = n_int
    dp = 2*pi/nphi
    A = np.array([0., 0., 0.])
    r = np.array([x, y, z])
    for ip in range(nphi):
        p = ip * dp
        rp = np.array([cos(p), sin(p), 0.])
        ep = np.array([-sin(p), cos(p), 0.])
        dist = np.linalg.norm(r - rp)
        tr = t - dist/c
        A += (tr>0) * sin(w*tr)/ dist * ep * dp
    return A / c**2 * a * I

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
    E = - diff_A_t(x, y, z, t, d)
    return E

########## B - magnetic field ##########

def B_compute(x, y, z, t, d):
    B = np.zeros(3)

    B[0] = diff_A_y(x, y, z, t, d)[2] - diff_A_z(x, y, z, t, d)[1]
    B[1] = diff_A_z(x, y, z, t, d)[0] - diff_A_x(x, y, z, t, d)[2]
    B[2] = diff_A_x(x, y, z, t, d)[1] - diff_A_y(x, y, z, t, d)[0]

    return B

########## Plot ##########

def plot(x, y, z, t0):
    t = np.linspace(t0, t0 + 8 * pi / w, 100)

    grid1 = np.zeros((len(t), 3))
    grid2 = np.zeros((len(t), 3))

    for i in range(len(t)):
        grid1[i,:] = E_compute(x, y, z, t[i], delta)
        grid2[i,:] = B_compute(x, y, z, t[i], delta)

    fig, ax1 = plt.subplots()
    ax1.plot(t, grid1[:, 0], label='$E_x$', color=(255/255,   0/255, 0/255), linestyle='solid')
    ax1.plot(t, grid1[:, 1], label='$E_y$', color=(255/255, 128/255, 0/255), linestyle='dashed')
    ax1.plot(t, grid1[:, 2], label='$E_z$', color=(255/255, 255/255, 0/255), linestyle='dotted')

    ax1.set_ylabel('E, <add units>')

    ax2 = ax1.twinx()
    ax2.plot(t, grid2[:, 0], label='$B_x$', color=(0/255,   0/255, 255/255), linestyle='solid')
    ax2.plot(t, grid2[:, 1], label='$B_y$', color=(0/255, 128/255, 255/255), linestyle='dashed')
    ax2.plot(t, grid2[:, 2], label='$B_z$', color=(0/255, 255/255, 255/255), linestyle='dotted')
    ax2.set_ylabel('B, <add units>')

    ax1.set_xlabel('t, s')
    fig.legend(loc='upper right')  # Adjust the legend position

    plt.title(f'Position: r = ({x}, {y}, {z}), $t_0$ = {t[0]}')

    # formatter = ticker.ScalarFormatter(useMathText=True)
    # formatter.set_scientific(True)
    # formatter.set_powerlimits((-1, 1))
    # ax1.yaxis.set_major_formatter(formatter)
    # ax2.yaxis.set_major_formatter(formatter)

    plt.tight_layout()
    plt.show()

########## Helper functions ##########

def distance(x, y, z, phi):
    return np.sqrt( (x - cos(phi))**2 + (y - sin(phi))**2 + z**2 )

########## Main ##########

def main():
    x, y, z, t0 = 10*a, 0, 0, 0
    plot(x, y, z, t0)
    T = np.linspace(0, 2*pi/w, 10)
    for t in T:
        A_compute(x, y, z, t)

##########################

if __name__ == '__main__':
    main()