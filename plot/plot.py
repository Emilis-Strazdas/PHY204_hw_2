import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
from matplotlib import ticker
from fields.e_field import E_compute
from fields.b_field import B_compute
from utils.initial_conditions import *


formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-1, 1))

def plot(x, y, z, t0):
    t = np.linspace(t0, t0 + 8 * pi / W, 100)

    grid1 = np.zeros((len(t), 3))
    grid2 = np.zeros((len(t), 3))

    for i in range(len(t)):
        grid1[i,:] = B_compute(x, y, z, t[i])
        grid2[i,:] = E_compute(x, y, z, t[i])

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('$t$, $s$')

    ax1.plot(t, grid1[:, 0], label='$B_x$', color=(0/255,   0/255, 255/255), linestyle='solid')
    ax1.plot(t, grid1[:, 1], label='$B_y$', color=(0/255, 128/255, 255/255), linestyle='dashed')
    ax1.plot(t, grid1[:, 2], label='$B_z$', color=(0/255, 255/255, 255/255), linestyle='dotted')
    ax1.set_ylabel('$B$, $G$')
    # ax1.yaxis.set_major_formatter(formatter)

    ax2 = ax1.twinx()
    ax2.plot(t, grid2[:, 0], label='$E_x$', color=(255/255,   0/255, 0/255), linestyle='solid')
    ax2.plot(t, grid2[:, 1], label='$E_y$', color=(255/255, 128/255, 0/255), linestyle='dashed')
    ax2.plot(t, grid2[:, 2], label='$E_z$', color=(255/255, 255/255, 0/255), linestyle='dotted')
    ax2.set_ylabel(r'$E$, $\frac{statV}{cm}$')
    # ax2.yaxis.set_major_formatter(formatter)

    fig.legend()  # Adjust the legend position
    plt.title(f'Position: r = ({x}, {y}, {z}), $t_0$ = {t[0]}')

    plt.tight_layout()
    plt.show()
