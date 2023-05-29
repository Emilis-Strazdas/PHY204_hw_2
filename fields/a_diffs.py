from fields.a_field import A_compute
from utils.calculus_constants import *


def diff_A_x(x, y, z, t):
    return (A_compute(x + d, y, z, t) - A_compute(x - d, y, z, t)) / (2 * d)

def diff_A_y(x, y, z, t):
    return (A_compute(x, y + d, z, t) - A_compute(x, y - d, z, t)) / (2 * d)

def diff_A_z(x, y, z, t):
    return (A_compute(x, y, z + d, t) - A_compute(x, y, z - d, t)) / (2 * d)

def diff_A_t(x, y, z, t):
    return (A_compute(x, y, z, t + d) - A_compute(x, y, z, t - d)) / (2 * d)