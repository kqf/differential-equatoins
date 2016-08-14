#!/usr/bin/python


import seaborn as sn
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

from Solver import Solver

# lorentz equations in terms of finite differences


def lorentz(uv, dx, dt):
    xn, yn, zn = uv
    sigma, rho, beta = 10, 28, 8. / 3.

    x = xn + dt * (yn - xn) * sigma
    y = yn + dt * (xn * (rho - zn) - yn)
    z = zn + dt * (xn * yn - beta * zn)
    return x, y, z


def get_lorentz(initial_conditions=(1., 0., 0.), step=0.01):
    solver = Solver((0, 2.), 0, 25, step)
    solver.initial_conditions = lambda x = 1: initial_conditions

    time = np.arange(0, 25, step)
    ans = solver.solve1d(lorentz)
    x, y, z = np.array(zip(*ans))
    return time, x, y, z


def main():
    fig = plt.figure(figsize=(15, 6))
    ax1 = fig.add_subplot(121, projection='3d')
    time, x, y, z = get_lorentz()
    ax1.plot(x, y, z, alpha=0.5)
    time, x, y, z = get_lorentz((-10., -10., -10.))
    ax1.plot(x, y, z, alpha=0.5)

    ax2 = fig.add_subplot(122)
    ax2.plot(time, x)
    ax2.plot(time, y)
    ax2.plot(time, z)
    plt.show()


if __name__ == '__main__':
    main()
