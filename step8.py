#!/usr/bin/python


import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

from Solver import Solver, MultiSolver


def profile(solver):
    u = np.ones(solver.nx)
    v = np.ones(solver.nx)
    dx, dy = solver.dx
    # set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
    u[.5 / dy:1 / dy + 1, .5 / dx:1 / dx + 1] = 2
    # set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
    v[.5 / dy:1 / dy + 1, .5 / dx:1 / dx + 1] = 2
    return u, v

# 2D convection


def convection(m, dr, dt):
    u, v = m
    un = u.copy()
    vn = v.copy()

    dx, dy = dr
    c = 1
    sigma = .0009
    nu = 0.01
    dt = sigma * dx * dy / nu

    u[1:-1, 1:-1] = un[1:-1, 1:-1] - dt / dx * un[1:-1, 1:-1] * (un[1:-1, 1:-1] - un[1:-1, 0:-2]) - dt / dy * vn[1:-1, 1:-1] * \
        (un[1:-1, 1:-1] - un[0:-2, 1:-1]) + nu * dt / dx**2 * (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) + \
        nu * dt / dy**2 * (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])

    v[1:-1, 1:-1] = vn[1:-1, 1:-1] - dt / dx * un[1:-1, 1:-1] * (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) - dt / dy * vn[1:-1, 1:-1] * \
        (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) + nu * dt / dx**2 * (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) + \
        nu * dt / dy**2 * (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])

    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1

    v[0, :] = 1
    v[-1, :] = 1
    v[:, 0] = 1
    v[:, -1] = 1

    return u, v


def main():
    msolver = MultiSolver((Solver((0, 2.), 41), Solver((0, 2.), 41)), 120, 0.1)
    msolver.initial_conditions = lambda x = None: profile(msolver)

    uinit, vinit = msolver.initial_conditions(msolver.argument())
    u, v = msolver.solve(convection)

    # the figsize parameter can be used to produce different sized images
    fig = plt.figure(figsize=(11, 7), dpi=100)
    ax = fig.gca(projection='3d')
    X, Y = np.meshgrid(*msolver.argument())

    ax.plot_wireframe(X, Y, u, cmap=cm.coolwarm)
    ax.plot_wireframe(X, Y, 2 - v, cmap=cm.coolwarm)
    plt.show()


if __name__ == '__main__':
    main()
