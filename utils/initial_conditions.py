"""
This file contains the initial conditions for the simulation.
"""

########## Initial conditions ##########

A = 1               # radius of the ring
I = 1               # current
W = 100             # angular velocity of the current
C = 10000           # speed of light

########## Point where the fields are calculated ##########

X, Y, Z, T0 = round(( 1 + 9 / 10000 * C)*A, 4), 0, 0, 0