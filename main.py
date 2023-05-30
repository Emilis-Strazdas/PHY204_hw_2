"""
Main file for the project.

Plots the fields at a given point in space and time.
"""

########## Imports ##########
from plot.plot import plot
from utils.initial_conditions import *

########## Main ##########

def main():
    plot(X, Y, Z, T0)

##############################

if __name__ == '__main__':
    main()