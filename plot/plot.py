"""
Plots the magnetic field B and the elctric field E at a point (x, y, z) as a
function of time t. The time range is 4 periods of the wave and starts at t0.
"""

########## Imports ##########

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
from matplotlib import ticker
from fields.e_field import E_compute
from fields.b_field import B_compute
from utils.initial_conditions import *

########## Functions ##########

# Set the scientific notation for the y-axes
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-1, 1))

def plot(x, y, z, t0):
    """
    Plots the magnetic field B and the elctric field E at a point (x, y, z) as a
    function of time t. The time range is 4 periods of the wave and starts at t0.

    Args:
        x, y, z: coordinates of the point where the fields are calculated
        t0: starting time of the plot

    Returns:
        None
    """
    t = np.linspace(t0, t0 + 5 * 2 * pi / W, 200)

    grid1 = np.zeros((len(t), 3))
    grid2 = np.zeros((len(t), 3))

    for i in range(len(t)):
        grid1[i,:] = B_compute(x, y, z, t[i])
        grid2[i,:] = E_compute(x, y, z, t[i])

    for i in range(len(t)):
        if np.linalg.norm(grid1[i,:]) > 0:
            grid1[i,:] = 0
            break

    for i in range(len(t)):
        if np.linalg.norm(grid2[i,:]) > 0:
            grid2[i,:] = 0
            break

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
    plt.title(f'Position: r = ({x}, {y}, {z}), $t_0$ = {t[0]}, c = {C}')

    plt.tight_layout()
    plt.show()
